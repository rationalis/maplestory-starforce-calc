use crate::consts::*;
use crate::prio::Prio;

use std::cmp::Ordering;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Mul, MulAssign};

use arrayfire::*;
use arrayfire::MatProp::NONE;

use rustc_hash::FxHashMap;

pub fn merge_or_insert<K, V, D>(dist: &mut FxHashMap<K, V>, key: K, p: D)
where
    K: Eq + Hash,
    V: AddAssign<D>,
    D: Copy + Into<V>,
{
    dist.entry(key)
        .and_modify(|p0| *p0 += p)
        .or_insert(p.into());
}

pub fn round_bucket(mesos: Meso) -> (u16, Meso) {
    if mesos <= IDENT_BINS as i32 {
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

#[derive(Clone)]
pub struct Distr {
    pub costs: Array<u32>,
    pub probs: Array<f64>,
}

impl Distr {
    fn new(costs: Array<u32>, probs: Array<f64>) -> Self {
        let mut res = Self { costs, probs };
        res.truncate(None);
        (res.costs, res.probs) = sort_by_key(&res.costs, &res.probs, 0, true);
        res
    }

    pub fn normalize(&mut self) {
        let total_prob: f64 = self.dist.iter().map(|(_, p)| *p).sum();
        self.dist.iter_mut().for_each(|(_, p)| *p /= total_prob);
    }

    pub fn truncate(&mut self, threshold: Option<f64>) {
        let threshold = match threshold {
            Some(p) => p,
            None => DIST_THRESHOLD,
        };
        (self.probs, self.costs) = sort_by_key(&self.probs, &self.costs, 0, false);
        // self.dist.sort_unstable_by_key(|(_, p)| f(*p));
        // self.dist.reverse();
        let mut total_prob = 0.0;
        let mut last_key = 0;
        for i in 0..self.dist.len() {
            total_prob += self.dist[i].1;
            if total_prob > 1.0 - threshold {
                last_key = self.dist[i].0;
                self.dist.truncate(i + 1);
                break;
            }
        }
        if total_prob < 1.0 {
            if last_key < 66 {
                self.dist.push((last_key, 1.0 - total_prob));
            } else {
                self.normalize();
            }
        }
    }

    pub fn constant(c: usize) -> Self {
        Self {
            dist: vec![(round_bucket(c as i32).0, 1.0)],
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
            dist.push((i, remaining * p));
            i += 1;
            remaining -= remaining * p;
        }
        Self { dist }
    }

    /// Calculate the joint distribution (i.e. negative multinomial) of outcomes.
    /// At checkpoint stars with 0 down chance (10, 15, 20), the number of downs
    /// is fixed to 0 and the number of attempts can vary, so this function will
    /// also map (downs, booms) to a distribution.
    pub fn sim(probs: [f64; 4]) -> FxHashMap<(u8, u8), PartialDistr> {
        let [up, stay, down, boom] = probs;
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
                merge_or_insert(&mut successes, (downs, booms), (attempts, p * up));
                update(&mut states, (downs + 1, booms), (attempts, p * down));
                update(&mut states, (downs, booms + 1), (attempts, p * boom));
                update(&mut states, (downs, booms), (attempts, p * stay));
            }
        }
        successes
    }

    pub fn shift(&mut self, c: Meso) -> &mut Self {
        self.dist
            .iter_mut()
            .for_each(|(c0, _)| *c0 = round_bucket(unbin(*c0) + c).0);
        self
    }

    pub fn scale(&mut self, n: Meso) -> &mut Self {
        self.dist
            .iter_mut()
            .for_each(|(c0, _)| *c0 = round_bucket(unbin(*c0) * n).0);
        self
    }

    pub fn add(&self, other: &Self) -> Self {
        if self.dist.len() == 1 && other.dist.len() == 1 {
            return Self {
                dist: vec![(BIN_SUMS[self.dist[0].0 as usize][other.dist[0].0 as usize], 1.0)]
            };
        }

        let (dist1_bin, dist1_prob): (Vec<_>, Vec<f64>) = self.dist.iter().cloned().unzip();
        let (dist2_bin, dist2_prob): (Vec<_>, Vec<f64>) = other.dist.iter().cloned().unzip();
        let dist1_bin = vec_to_af(&dist1_bin);
        let dist2_bin = vec_to_af(&dist2_bin);
        let dist1_prob = vec_to_af(&dist1_prob);
        let dist2_prob = transpose(&vec_to_af(&dist2_prob), false);
        let prod = matmul(&dist1_prob, &dist2_prob, NONE, NONE);
        let bin_sums_mat = moddims(&BIN_SUMS_AF,
                                   dim4!(NUM_BINS as u64,
                                         NUM_BINS as u64, 1, 1));
        let out_bins = view!(bin_sums_mat[dist1_bin, dist2_bin]);
        let keys = flat(&out_bins);
        let prod = flat(&prod);
        let (keys, prod) = if prod.elements() > 1 {
            let sorted = sort_by_key(&keys, &prod, 0, true);
            // let sorted = (view!(BIN_SUMS_AF[PERMUTE.0]), view!(prod[PERMUTE.1]));
            let sorted_keys: Array<u32> = sorted.0.cast();
            let reduced = sum_by_key(&sorted_keys, &sorted.1, 0);
            reduced
        } else {
            (keys.cast(), prod)
        };

        let dist_af =
            af_to_vec(&keys).into_iter()
            .zip(af_to_vec(&prod))
            .map(|(i, p)| (i as u16, p))
            .collect();
        Self::new(dist_af)
    }

    /// Calculate the mean and standard deviation.
    pub fn stats(&self) -> (u64, u64) {
        let mean: f64 = self.dist.iter().map(|(c, p)| (unbin(*c) as f64) * p).sum();
        let stddev: f64 = self
            .dist
            .iter()
            .map(|(c, p)| (unbin(*c) as f64 - mean).powi(2) * p)
            .sum::<f64>()
            .sqrt();
        let mean = (UNIT as f64) * mean;
        let stddev = (UNIT as f64) * stddev;
        (mean as u64, stddev as u64)
    }

    /// Calculate the lowest known cost, 3 quartiles, and the 99th percentile.
    pub fn quartiles(&self) -> (u64, u64, u64, u64, u64) {
        let mut cdf = 0.0;
        let dist = &self.dist;
        let mut i = 0;
        let mut quartile1 = 0;
        while cdf <= 0.25 {
            quartile1 = dist[i].0;
            cdf += dist[i].1;
            i += 1;
        }
        let mut quartile2 = 0;
        while cdf <= 0.5 {
            quartile2 = dist[i].0;
            cdf += dist[i].1;
            i += 1;
        }
        let mut quartile3 = 0;
        while cdf <= 0.75 {
            quartile3 = dist[i].0;
            cdf += dist[i].1;
            i += 1;
        }
        let mut max = 0;
        while cdf <= 0.99 {
            max = dist[i].0;
            cdf += dist[i].1;
            i += 1;
        }
        let min = dist.first().unwrap().0;

        (
            unbin(min) as u64 * UNIT as u64,
            unbin(quartile1) as u64 * UNIT as u64,
            unbin(quartile2) as u64 * UNIT as u64,
            unbin(quartile3) as u64 * UNIT as u64,
            unbin(max) as u64 * UNIT as u64,
        )
    }

    pub fn iter(&self) -> impl Iterator<Item = (Meso, f64)> + '_ {
        self.dist.iter().map(|&(i, p)| (unbin(i), p))
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
                merge_or_insert(&mut res.dist, round_bucket(c + c2).1, p * p2);
                res.total += p * p2;
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
        Distr::new(
            self.dist
                .into_iter()
                .map(|(c, p)| (round_bucket(c).0, g(p)))
                .collect(),
        )
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
        res.dist = self.dist.iter().map(|(k, v)| (k * other, *v)).collect();
        res.total = self.total;
        res
    }
}
