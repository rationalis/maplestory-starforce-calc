use crate::consts::*;
use crate::prio::Prio;

use std::cmp::Ordering;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Mul, MulAssign};

use rustc_hash::FxHashMap;

#[derive(Clone, Debug)]
pub struct Distr {
    pub dist: Vec<(Meso, f64)>
}

pub fn slice_to_distr(slice: &[f64]) -> Distr {
    Distr {
        dist: slice
            .into_iter()
            .enumerate()
            .map(|(n, p)| (n as i32, *p))
            .collect()
    }
}

pub fn merge_or_insert<K, V, D>(dist: &mut FxHashMap<K, V>, key: K, p: D) where
    K: Eq + Hash, V: AddAssign<D>, D: Copy + Into<V>
{
    dist.entry(key)
        .and_modify(|p0| *p0 += p)
        .or_insert(p.into());
}

pub fn round_bucket(mesos: Meso) -> Meso {
    if mesos <= 1000 {
        return mesos;
    }
    let c = mesos.abs();
    let res = BINS.binary_search(&c);
    let idx = match res {
        Ok(found) => found,
        Err(insertion) => {
            let mut closest = None;
            let mut closest_diff = None;
            if insertion < BINS.len() {
                closest = Some(insertion);
                closest_diff = Some((BINS[insertion] - c).abs());
            }
            if insertion > 0 {
                if closest.is_none() || (BINS[insertion - 1] - c).abs() < closest_diff.unwrap() {
                    closest = Some(insertion - 1);
                }
            }
            closest.unwrap()
        }
    };
    let c = BINS[idx];

    c * mesos.signum()
}

impl Distr {
    pub fn new(dist: Vec<(Meso, f64)>) -> Self {
        let mut map: FxHashMap<_, _> = Default::default();
        for (c, p) in dist.into_iter() {
            merge_or_insert(&mut map, round_bucket(c), p);
        }
        let dist: Vec<_> = map.into_iter().collect();
        let mut res = Self { dist };
        res.truncate(None);
        // if res.dist.len() > 14000 {
        //     dbg!("Distr longer than 14000");
        //     dbg!(&res.dist.iter().map(|(p,_)| *p).collect::<Vec<_>>()[0..1000]);
        //     dbg!(&res.dist[0..100]);
        //     for &(c, p) in res.dist.iter() {
        //         assert_eq!(c, round_bucket(c));
        //         if c >= 100_000_000 {
        //             println!("{},{}", c, p);
        //         }
        //     }
        //     panic!();
        // }
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
        let mut last_key = 0;
        for i in 0..self.dist.len() {
            total_prob += self.dist[i].1;
            if total_prob > 1.0 - threshold {
                last_key = self.dist[i].0;
                self.dist.truncate(i+1);
                break;
            }
        }
        if total_prob < 1.0 {
            if last_key < 66 {
                self.dist.push((last_key, 1.0 - total_prob));
            } else {
                // self.dist.iter_mut().for_each(|(_,p)| *p /= total_prob);
                self.normalize();
            }
        }
    }

    pub fn constant(c: Meso) -> Self {
        Self {
            dist: vec![(c, 1.0)],
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

    pub fn chance_time(&self, normal_up_cost: &Self, chance_time_cost: i32) -> Self {
        let mut dist = FxHashMap::default();
        for &(n, p) in self.dist.iter() {
            if n == 0 {
                dist.insert(0, p);
                continue;
            }
            if n == 1 {
                dist.insert(chance_time_cost, p);
                continue;
            }
            for (c, p2) in normal_up_cost.dist.iter() {
                merge_or_insert(&mut dist, round_bucket((n-1)*c + chance_time_cost), p*p2);
            }
        }
        let dist = dist.into_iter().collect();
        dbg!(self.dist.len(), normal_up_cost.dist.len());
        Self::new(dist)
    }

    /// This method is actually faster by something like 25x compared to summing
    /// the probabilities over the compound binomial distribution.
    pub fn booms(succ_rate: f64, boom_rate: f64) -> Self {
        let mut states: [f64; MAX_BOOMS] = [0.0; MAX_BOOMS];
        let mut successes: [f64; MAX_BOOMS] = [0.0; MAX_BOOMS];
        states[0] = 1.0;
        let mut total_prob = 0.0;
        while total_prob < 1.0 - DIST_THRESHOLD {
            debug_assert!((total_prob + states.iter().sum::<f64>() - 1.0).abs() < 1e-6);
            for booms in 0..MAX_BOOMS {
                let p = states[booms];
                if p > 1e-14 {
                    successes[booms] += p * succ_rate;
                    total_prob += p * succ_rate;
                    states[booms] = p * (1.0 - succ_rate - boom_rate);
                    states[booms+1] += p * boom_rate;
                    break;
                }
            }
        }
        slice_to_distr(&successes)
    }

    pub fn shift(&mut self, c: Meso) -> &mut Self {
        self.dist.iter_mut().for_each(|(c0, _)| *c0 = round_bucket(*c0 + c));
        self
    }

    pub fn scale(&mut self, n: Meso) -> &mut Self {
        self.dist.iter_mut().for_each(|(c0, _)| *c0 = round_bucket(*c0 * n));
        self
    }

    pub fn product(&self, other: &Self) -> Self {
        let mut dist = FxHashMap::default();
        for (i, p) in self.dist.iter() {
            for (c, p2) in other.dist.iter() {
                merge_or_insert(&mut dist, round_bucket(i*c), p*p2);
            }
        }
        let dist = dist.into_iter().collect();
        Self::new(dist)
    }

    pub fn add(&self, other: &Self) -> Self {
        let mut dist = FxHashMap::default();
        for (c, p) in self.dist.iter() {
            for (c2, p2) in other.dist.iter() {
                merge_or_insert(&mut dist, round_bucket(c+c2), p*p2);
            }
        }
        let dist: Vec<_> = dist.into_iter().collect();
        Self::new(dist)
    }

    pub fn mix(&mut self, other: &PartialDistr) -> &mut Self {
        panic!()
    }

    pub fn expected_cost(&self) -> u64 {
        let sum: f64 = self.dist.iter().map(|(c, p)| (*c as f64)*p).sum();
        let sum = (UNIT as f64) * sum;
        sum as u64
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
        *self = self.add(other)
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
        let res = self.clone();
        res.scale(other);
        res
    }
}

impl Add<PartialDistr> for Distr {
    // stuff
}

#[derive(Debug, Default, Eq, PartialEq)]
pub struct PartialDistr {
    total: F,
    dist: FxHashMap<Meso, F>,
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
        let mut res = *self.clone();
        res.dist.iter_mut().for_each(|(c, _)| *c *= other);
        res
    }
}

