use crate::binning::*;
use crate::consts::*;
use crate::prio::Prio;

use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Mul, MulAssign};

use rustc_hash::FxHashMap;

#[derive(Clone, Debug)]
pub struct Distr {
    pub dist: Vec<(Bin, f64)>,
}

impl Distr {
    fn new(dist: Vec<(Bin, f64)>) -> Self {
        let mut res = Self { dist };
        res.truncate(None);
        res.dist.sort_unstable_by_key(|(c, _)| *c);

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
        self.dist.sort_unstable_by_key(|(_, p)| f(*p));
        // if there is any dupe, dedupe
        for i in 1..self.dist.len() {
            if self.dist[i-1].0.raw() == self.dist[i].0.raw() {
                let mut map = FxHashMap::default();
                for &(c, p) in self.dist.iter() {
                    merge_or_insert(&mut map, c, p);
                }
                self.dist = map.into_iter().collect();
                self.dist.sort_unstable_by_key(|(_, p)| f(*p));
                break;
            }
        }
        self.dist.reverse();
        let mut total_prob = 0.0;
        let mut last_key = 0;
        for i in 0..self.dist.len() {
            total_prob += self.dist[i].1;
            if total_prob > 1.0 - threshold {
                last_key = self.dist[i].0.raw();
                self.dist.truncate(i + 1);
                break;
            }
        }
        if total_prob < 1.0 {
            if last_key < 66 {
                self.dist.push((Bin::from_id(last_key), 1.0 - total_prob));
            } else {
                self.normalize();
            }
        }
    }

    pub fn constant(c: Meso) -> Self {
        Self {
            dist: vec![(Bin::from_cost(c), 1.0)],
        }
    }

    pub fn zero() -> Self {
        Self::constant(f(0))
    }

    pub fn geom(p: f64) -> Self {
        let mut dist = Vec::new();
        let mut remaining = 1.0;
        let mut i = 1;
        while remaining > DIST_THRESHOLD {
            dist.push((Bin::from_int(i), remaining * p));
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
        states.push((0, 0), (f(0), f(1)));
        while states.total_prob > 1e-6 {
            let ((downs, booms), pdist) = states.pop();
            for (attempts, &p) in pdist.dist.iter() {
                let attempts = *attempts + 1.0;
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
            .for_each(|(c0, _)| *c0 = Bin::from_cost(c0.unbin() + c));
        self
    }

    pub fn scale(&mut self, n: Meso) -> &mut Self {
        self.dist
            .iter_mut()
            .for_each(|(c0, _)| *c0 = Bin::from_cost(c0.unbin() * n));
        self
    }

    pub fn add(&self, other: &Self) -> Self {
        // let mut dist = [0.0; NUM_BINS];
        let mut dist = vec![0.0; ADAPTIVE_BINS.len()];
        for &(c, p) in self.dist.iter() {
            // let lookup = &BIN_SUMS[c.raw() as usize];
            let row_offset = c.raw() as usize * ADAPTIVE_BINS.len();
            for &(c2, p2) in other.dist.iter() {
                // let i = lookup[c2.raw() as usize];
                let i = BIN_SUMS[row_offset + c2.raw() as usize];
                dist[i as usize] += p * p2;
            }
        }
        let dist = dist
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| if p == 0.0 { None } else { Some((Bin::from_id(i as BinId), p)) })
            .collect();
        Self::new(dist)
    }

    /// Calculate the mean and standard deviation.
    pub fn stats(&self) -> (Meso, Meso) {
        let mean: Meso = self.dist.iter().map(|(c, p)| (c.unbin()) * p).sum();
        let stddev = self
            .dist
            .iter()
            .map(|(c, p)| g(c.unbin() - mean).powi(2) * p)
            .sum::<f64>()
            .sqrt();
        // let mean = (UNIT as f64) * mean;
        // let stddev = (UNIT as f64) * stddev;
        (mean, f(stddev))
    }

    /// Calculate the lowest known cost, 3 quartiles, and the 99th percentile.
    pub fn quartiles(&self) -> (Meso, Meso, Meso, Meso, Meso) {
        let mut cdf = 0.0;
        let dist = &self.dist;
        let mut i = 0;
        let mut quartile1 = f(0);
        while cdf <= 0.25 {
            quartile1 = dist[i].0.unbin();
            cdf += dist[i].1;
            i += 1;
        }
        let mut quartile2 = f(0);
        while cdf <= 0.5 {
            quartile2 = dist[i].0.unbin();
            cdf += dist[i].1;
            i += 1;
        }
        let mut quartile3 = f(0);
        while cdf <= 0.75 {
            quartile3 = dist[i].0.unbin();
            cdf += dist[i].1;
            i += 1;
        }
        let mut max = f(0);
        while cdf <= 0.99 {
            max = dist[i].0.unbin();
            cdf += dist[i].1;
            i += 1;
        }
        let min = dist.first().unwrap().0.unbin();

        (
            min, quartile1, quartile2, quartile3, max,
            // unbin(min) as u64 * UNIT as u64,
            // unbin(quartile1) as u64 * UNIT as u64,
            // unbin(quartile2) as u64 * UNIT as u64,
            // unbin(quartile3) as u64 * UNIT as u64,
            // unbin(max) as u64 * UNIT as u64,
        )
    }

    pub fn iter(&self) -> impl Iterator<Item = (Meso, f64)> + '_ {
        self.dist.iter().map(|&(i, p)| (i.unbin(), p))
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

impl Add<Meso> for &Distr {
    type Output = Distr;

    fn add(self, other: Meso) -> Self::Output {
        let mut res = self.clone();
        res.shift(other);
        res
    }
}

impl MulAssign<Meso> for Distr {
    fn mul_assign(&mut self, other: Meso) {
        self.scale(other);
    }
}

impl Mul<Meso> for &Distr {
    type Output = Distr;

    fn mul(self, other: Meso) -> Self::Output {
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
                // merge_or_insert(&mut res.dist, round_bucket(c + c2).1, p * p2);
                merge_or_insert(&mut res.dist, c + c2, p * p2);
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
                .map(|(c, p)| (Bin::from_cost(c), g(p)))
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

impl Mul<Meso> for &PartialDistr {
    type Output = PartialDistr;

    fn mul(self, other: Meso) -> Self::Output {
        let mut res = PartialDistr::default();
        res.dist = self.dist.iter().map(|(k, v)| (*k * other, *v)).collect();
        res.total = self.total;
        res
    }
}
