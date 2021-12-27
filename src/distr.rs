use crate::consts::*;
use crate::prio::Prio;

use std::cmp::{Ordering, PartialEq};
use std::hash::Hash;
use std::ops::{Add, AddAssign, Mul, MulAssign};

use ndarray::{Array, Array1};
use rustc_hash::FxHashMap;

pub fn round_bucket(mesos: Meso) -> (u16, Meso) {
    let repr: f32 = mesos.into();
    let idx = repr.log(BIN_EXP as f32) - BASE_STAR_FACTOR.log(BIN_EXP as f32);
    let idx = idx.round() as u16;
    (idx, BINS[idx as usize])
    //round_bucket_impl(&*BINS, mesos)
}

pub fn unbin(u: u16) -> Meso {
    BINS[u as usize]
}

pub type Distr = Array1<f64>;

pub fn zero_distr() -> Distr {
    Array::zeros((NUM_BINS))
}

pub fn geom() -> Distr {
    let mut dist = Vec::new();
    let mut remaining = 1.0;
    let mut i = 1;
    while remaining > DIST_THRESHOLD {
        dist.push((i, remaining * p));
        i += 1;
        remaining -= remaining * p;
    }
    Array::from_vec(dist)
}

pub fn fftlog(d: Distr) -> Distr {
}

/*
#[derive(Clone, Debug)]
pub struct Distr {
    pub dist: Vec<(u16, f64)>,
}

impl Distr {
    fn new(dist: Vec<(u16, f64)>) -> Self {
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
        self.dist.reverse();
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
    pub fn sim(probs: [f64; 4]) -> FxHashMap<(u8, u8), PartialDistr<i32>> {
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
        let mut dist = [0.0; NUM_BINS];
        for &(c, p) in self.dist.iter() {
            let lookup = &BIN_SUMS[c as usize];
            for &(c2, p2) in other.dist.iter() {
                let i = lookup[c2 as usize];
                dist[i as usize] += p * p2;
            }
        }
        let dist = dist
            .iter()
            .enumerate()
            .filter_map(|(i, &p)| if p == 0.0 { None } else { Some((i as u16, p)) })
            .collect();
        Self::new(dist)
    }

    /// Calculate the mean and standard deviation.
    pub fn stats(&self) -> (u64, u64) {
        let mean: f64 = self.dist.iter().map(|(c, p)| (unbin(*c).into()) * p).sum();
        let stddev: f64 = self
            .dist
            .iter()
            .map(|(c, p)| (f64::from(unbin(*c)) - mean).powi(2) * p)
            .sum::<f64>()
            .sqrt();
        // TODO
        // let mean = (UNIT as f64) * mean;
        // let stddev = (UNIT as f64) * stddev;
        // (mean as u64, stddev as u64)

        (0, 0)
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
            0,0,0,0,0
            // TODO
            // unbin(min) as u64 * UNIT as u64,
            // unbin(quartile1) as u64 * UNIT as u64,
            // unbin(quartile2) as u64 * UNIT as u64,
            // unbin(quartile3) as u64 * UNIT as u64,
            // unbin(max) as u64 * UNIT as u64,
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

impl Add<&PartialDistr<Meso>> for Distr {
    type Output = PartialDistr<Meso>;

    fn add(self, other: &PartialDistr<Meso>) -> Self::Output {
        let mut res = PartialDistr::default();
        for (c, p) in self.iter() {
            let p = f(p);
            for (c2, p2) in other.dist.iter() {
                merge_or_insert(&mut res.dist, round_bucket(c + *c2).1, p * p2);
                res.total += p * p2;
            }
        }
        res
    }
}

impl Into<Vec<(u64, f64)>> for Distr {
    fn into(self) -> Vec<(u64, f64)> {
        self.dist.iter()
            // .map(|&(c, p)| (unbin(c) as u64 * UNIT as u64, p))
            // TODO
            .map(|&(c, p)| (0u64, 0f64))
            .collect()
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct PartialDistr<T>
where T: Hash + PartialEq
{
    pub total: F,
    pub dist: FxHashMap<T, F>,
}

impl<T> PartialDistr<T> where T: Hash + PartialEq {
    pub fn mix(&mut self, other: Self) {
        for (&c, &p) in other.dist.iter() {
            merge_or_insert(&mut self.dist, c, p);
            self.total += p;
        }
    }
}

impl<T> Into<Distr> for PartialDistr<T> where T: Hash + PartialEq {
    fn into(self) -> Distr {
        Distr::new(
            self.dist
                .into_iter()
                .map(|(c, p)| (round_bucket(c).0, g(p)))
                .collect(),
        )
    }
}

impl<T> Ord for PartialDistr<T> where T: Hash + PartialEq + Ord
{
    fn cmp(&self, other: &Self) -> Ordering {
        self.total.cmp(&other.total)
    }
}

impl<T> PartialOrd for PartialDistr<T> where T: Hash + PartialEq + PartialOrd
{
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.total.cmp(&other.total))
    }
}

impl<T> From<(Meso, F)> for PartialDistr<T> where T: Hash + PartialEq {
    fn from(val: (Meso, F)) -> Self {
        let mut res = Self::default();
        res += val;
        res
    }
}

impl<T> From<&PartialDistr<T>> for F where T: Hash + PartialEq {
    fn from(distr: &PartialDistr<T>) -> F {
        distr.total
    }
}

impl<T> AddAssign<(T, F)> for PartialDistr<T>
    where T: PartialEq + Eq + Hash
{
    fn add_assign(&mut self, other: (T, F)) {
        let (c, p) = other;
        merge_or_insert(&mut self.dist, c, p);
        self.total += p;
    }
}

impl<T> Mul<T> for &PartialDistr<T>
    where T: PartialEq + Hash + Mul, <T as Mul>::Output: Default
{
    type Output = PartialDistr<T>;

    fn mul(self, other: T) -> Self::Output {
        let mut res = PartialDistr::default();
        res.dist = self.dist.iter().map(|(k, v)| (*k * other, *v)).collect();
        res.total = self.total;
        res
    }
}

*/