use crate::consts::*;
use crate::prio::Prio;

use std::cmp::Ordering;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Mul, MulAssign};

use rustc_hash::FxHashMap;

pub fn merge_or_insert<K, V, D>(dist: &mut FxHashMap<K, V>, key: K, p: D) where
    K: Eq + Hash, V: AddAssign<D>, D: Copy + Into<V>
{
    dist.entry(key)
        .and_modify(|p0| *p0 += p)
        .or_insert(p.into());
}

pub fn round_bucket(mesos: Meso) -> (u16, Meso) {
    if mesos <= 1000 {
        return (mesos as u16, mesos);
    }
    round_bucket_impl(&*BINS, mesos)
}

pub fn round_bucket_impl(bins: &[Meso], mesos: Meso) -> (u16, Meso) {
    let c = mesos;
    let res = bins.binary_search(&c);
    let idx = match res {
        Ok(found) => found,
        Err(insertion) => {
            let mut closest = None;
            let mut closest_diff = None;
            if insertion < bins.len() {
                closest = Some(insertion);
                closest_diff = Some((bins[insertion] - c).abs());
            }
            if insertion > 0 {
                if closest.is_none() || (bins[insertion - 1] - c).abs() < closest_diff.unwrap() {
                    closest = Some(insertion - 1);
                }
            }
            closest.unwrap()
        }
    };
    (idx as u16, bins[idx])
}

pub fn unbin(u: u16) -> Meso {
    BINS[u as usize]
}

#[derive(Clone, Debug)]
pub struct Distr {
    pub dist: Vec<(u16, f64)>
}

impl Distr {
    fn new(dist: Vec<(u16, f64)>) -> Self {
        let mut res = Self { dist };
        res.truncate(None);
        res
    }

    pub fn normalize(&mut self) {
        let total_prob: f64 = self.dist.iter().map(|(_, p)| *p).sum();
        self.dist.iter_mut().for_each(|(_,p)| *p /= total_prob);
    }

    pub fn truncate(&mut self, threshold: Option<f64>) {
        let threshold = match threshold {
            Some(p) => p,
            None => DIST_THRESHOLD
        };
        self.dist.sort_unstable_by_key(|(_,p)| f(*p));
        self.dist.reverse();
        let mut total_prob = 0.0;
        for i in 0..self.dist.len() {
            total_prob += self.dist[i].1;
            if total_prob > 1.0 - threshold {
                self.dist.truncate(i+1);
                break;
            }
        }
        self.dist.sort_unstable_by_key(|(c, _)| *c);
        let max = self.dist.last().unwrap().0;
        if total_prob < 1.0 {
            // self.dist.push((max + 1, 1.0 - total_prob));
            if max < 66 {
                self.dist.push((max+1, 1.0 - total_prob));
            } else {
                self.normalize();
            }
        }
    }

    pub fn constant(c: Meso) -> Self {
        Self {
            dist: vec![(round_bucket(c).0, 1.0)],
        }
    }

    pub fn zero() -> Self {
        Self::constant(0)
    }

    pub fn geom(p: f64) -> Self {
        let mut dist = Vec::new();
        let mut remaining = 1.0;
        let mut i = 1;
        while remaining > DIST_THRESHOLD {
            dist.push((i, remaining*p));
            i += 1;
            remaining -= remaining * p;
        }
        Self {
            dist
        }
    }

    /// Calculate the joint distribution (i.e. negative multinomial) of outcomes.
    /// At checkpoint stars with 0 down chance (10, 15, 20), the number of downs
    /// is fixed to 0 and the number of attempts can vary, so this function will
    /// also map (downs, booms) to a distribution.
    pub fn sim(start: Star) -> FxHashMap<(u8, u8), PartialDistr> {
        let [up, stay, down, boom] = PROBS[(start - 10) as usize];
        let mut successes = FxHashMap::default();
        let mut states: Prio<(u8, u8), PartialDistr> = Prio::new();
        let update = |states: &mut Prio<_, _>, k, v: (Meso, F)| {
            if v.1 < f(1e-16) {
                return;
            }
            states.push(k, v);
        };
        states.push((0, 0), (0, f(1.0)));
        while states.total_prob > 1e-6 {
            let ((downs, booms), pdist) = states.pop();
            for (attempts, &p) in pdist.dist.iter() {
                let attempts = attempts + 1;
                merge_or_insert(&mut successes, (downs, booms), (attempts, p*up));
                update(&mut states, (downs+1, booms), (attempts, p*down));
                update(&mut states, (downs, booms+1), (attempts, p*boom));
                update(&mut states, (downs, booms), (attempts, p*stay));
            }
        }
        successes
    }

    pub fn shift(&mut self, c: Meso) -> &mut Self {
        self.dist.iter_mut().for_each(|(c0, _)| *c0 = round_bucket(unbin(*c0) + c).0);
        self
    }

    pub fn scale(&mut self, n: Meso) -> &mut Self {
        self.dist.iter_mut().for_each(|(c0, _)| *c0 = round_bucket(unbin(*c0) * n).0);
        self
    }

    pub fn add(&self, other: &Self) -> Self {
        use ndarray::{Array1, Array2, Array, ArrayBase, Ix2};
        // type V = ArrayD<f64>;
        let dist1: (Vec<_>, Vec<f64>) = self.dist.iter().cloned().unzip();
        let dist1: Array1<_> = dist1.1.into();
        let dist1: Array2<_> = dist1.into_shape((self.dist.len(), 1)).unwrap();
        let dist2: (Vec<_>, Vec<f64>) = other.dist.iter().cloned().unzip();
        let dist2: Array1<_> = dist2.1.into();
        let dist2: Array2<_> = dist2.into_shape((1, other.dist.len())).unwrap();

        let outer_product = dist1.dot(&dist2);
        let mut dist = [0.0; NUM_BINS];
        for (i, row) in outer_product.outer_iter().enumerate() {
            let lookup = &BIN_SUMS[self.dist[i].0 as usize];
            for (j, elem) in row.indexed_iter() {
                let i = lookup[other.dist[j].0 as usize];
                dist[i as usize] += elem;
            }
        }

        // for &(c, p) in self.dist.iter() {
        //     let lookup = &BIN_SUMS[c as usize];
        //     let p = p;
        //     for &(c2, p2) in other.dist.iter() {
        //         let i = lookup[c2 as usize] as usize;
        //         dist[i] += p*p2;
        //     }
        // }
        let dist = dist
            .iter()
            .enumerate()
            .filter_map(
                |(i, &p)|
                {
                    if p == 0.0 {
                        None
                    } else {
                        Some((i as u16, p))
                    }
                }
            ).collect();
        Self::new(dist)
    }

    pub fn expected_cost(&self) -> u64 {
        let sum: f64 = self.dist.iter().map(|(c, p)| (unbin(*c) as f64)*p).sum();
        let sum = (UNIT as f64) * sum;
        sum as u64
    }

    pub fn iter(&self) -> impl Iterator<Item=(Meso, f64)> + '_ {
        self.dist.iter().map(|&(i, p)| (unbin(i), p))
    }
}

impl Into<Vec<(u64, f64)>> for Distr {
    fn into(self) -> Vec<(u64, f64)> {
        self.dist.iter()
            .map(|&(c, p)| (unbin(c) as u64 * UNIT as u64, p))
            .collect()
    }
}

impl Add for &Distr {
    type Output = Distr;

    fn add(self, other: Self) -> Self::Output {
        Distr::add(self, other)
    }
}

impl AddAssign<&Distr> for Distr {
    fn add_assign(&mut self, other: &Self) {
        *self = Distr::add(self, other)
    }
}

impl Add<i32> for &Distr {
    type Output = Distr;

    fn add(self, other: i32) -> Self::Output {
        let mut res = self.clone();
        res.shift(other);
        res
    }
}

impl MulAssign<i32> for Distr {
    fn mul_assign(&mut self, other: i32) {
        self.scale(other);
    }
}

impl Mul<i32> for &Distr {
    type Output = Distr;

    fn mul(self, other: i32) -> Self::Output {
        let mut res = self.clone();
        res.scale(other);
        res
    }
}

impl Add<&PartialDistr> for Distr {
    type Output = PartialDistr;

    fn add(self, other: &PartialDistr) -> Self::Output {
        let mut res = PartialDistr::default();
        for (c, p) in self.iter() {
            let p = f(p);
            for (c2, p2) in other.dist.iter() {
                merge_or_insert(&mut res.dist, round_bucket(c+c2).1, p*p2);
                res.total += p*p2;
            }
        }
        res
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct PartialDistr {
    pub total: F,
    pub dist: FxHashMap<Meso, F>,
}

impl PartialDistr {
    pub fn mix(&mut self, other: Self) {
        for (&c, &p) in other.dist.iter() {
            merge_or_insert(&mut self.dist, c, p);
            self.total += p;
        }
    }
}

impl Into<Distr> for PartialDistr {
    fn into(self) -> Distr {
        Distr::new(self.dist.into_iter().map(|(c,p)| (round_bucket(c).0,g(p))).collect())
    }
}

impl Ord for PartialDistr {
    fn cmp(&self, other: &Self) -> Ordering {
        self.total.cmp(&other.total)
    }
}

impl PartialOrd for PartialDistr {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.total.cmp(&other.total))
    }
}

impl From<(Meso, F)> for PartialDistr {
    fn from(val: (Meso, F)) -> Self {
        let mut res = Self::default();
        res += val;
        res
    }
}

impl From<&PartialDistr> for F {
    fn from(distr: &PartialDistr) -> F {
        distr.total
    }
}

impl AddAssign<(Meso, F)> for PartialDistr {
    fn add_assign(&mut self, other: (Meso, F)) {
        let (c, p) = other;
        merge_or_insert(&mut self.dist, c, p);
        self.total += p;
    }
}

impl Mul<i32> for &PartialDistr {
    type Output = PartialDistr;

    fn mul(self, other: i32) -> Self::Output {
        let mut res = PartialDistr::default();
        res.dist = self.dist.iter().map(|(k, v)| (k*other, *v)).collect();
        res.total= self.total;
        res
    }
}

