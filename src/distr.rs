use crate::consts::*;
use crate::prio::Prio;

use std::hash::Hash;

use rustc_hash::FxHashMap;

#[derive(Clone, Debug)]
pub struct Distr {
    pub dist: Vec<(Meso, f64)>
}

pub fn slice_to_distr(slice: &[f64]) -> Distr {
    Distr::new(slice.into_iter().enumerate().map(|(n, p)| (n as i32, *p)).collect())
}

pub fn merge_or_insert<K: Eq + Hash>(dist: &mut FxHashMap<K, f64>, key: K, p: f64) {
    dist.entry(key)
        .and_modify(|p0| *p0 += p)
        .or_insert(p);
}

pub fn round_bucket(mesos: Meso) -> Meso {
    let c = mesos;
    if c > 1_000_000 { // 100 bil
        round(c, 10_000)
    } else if c > 100_000 { // 10 bil
        round(c, 1000)
    } else if c > 10000 { // bil
        round(c, 100)
    } else if c > 1000 {
        round(c, 10)
    } else {
        c
    }
}

impl Distr {
    pub fn new(dist: Vec<(Meso, f64)>) -> Self {
        let mut map: FxHashMap<_, _> = Default::default();
        for (c, p) in dist.into_iter() {
            merge_or_insert(&mut map, round_bucket(c), p);
        }
        let dist: Vec<_> = map.into_iter().collect();
        let mut res = Self { dist };
        res.truncate();
        res
    }

    pub fn truncate(&mut self) {
        self.dist.sort_unstable_by_key(|(_,p)| f(*p));
        self.dist.reverse();
        let mut total_prob = 0.0;
        let mut last_key = 0;
        for i in 0..self.dist.len() {
            total_prob += self.dist[i].1;
            if total_prob > 1.0 - DIST_THRESHOLD {
                last_key = self.dist[i].0;
                self.dist.truncate(i+1);
                break;
            }
        }
        if total_prob < 1.0 {
            if last_key < 66 {
                self.dist.push((last_key, 1.0 - total_prob));
            } else {
                self.dist.iter_mut().for_each(|(_,p)| *p /= total_prob);
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

    pub fn downs(start: Star) -> [[f64; MAX_DOWNS]; 4] {
        type State = (u8, u8, u8, u8, bool);
        let [up, stay, down, boom] = PROBS_F64[(start - 10) as usize];
        let [downup, downstay, downdown, downboom] = PROBS_F64[(start - 11) as usize];
        let mut states: Prio<State, (F, i32)> = Prio::new();
        states.prob = Some(|v: &mut (F, i32)| {
            let out: &mut F = &mut v.0;
            out
        });
        let update = |map: &mut Prio<State, (F, i32)>, elem: (F, u8, u8, u8, u8, bool)| {
            let (p, du, ds, dd, db, at_start) = elem;
            if p == f(0.0) {
                return;
            }
            map.push((du, ds, dd, db, at_start), (p, -((du + ds + dd + db) as i32)));
        };
        states.push((0,0,0,0,true), (f(1.0), 0));
        let mut output: [[f64; MAX_DOWNS]; 4] = [[0.0; MAX_DOWNS]; 4];
        let mut joint: FxHashMap<(u8, u8, u8, u8), f64> = Default::default();
        let mut total_prob = 0.0;
        let mut uniq_subdists = 0;
        let mut distr_adds = 0;
        while total_prob < 1.0 - DIST_THRESHOLD {
            let ((du, ds, dd, db, at_start), (p, _)) = states.pop();
            if p < 1e-12 {
                continue;
            }
            let states = &mut states;
            if at_start {
                let mut remaining = g(p);
                let mut succ = 0.0;
                let mut p_down = f(0.0);
                while remaining > 1e-16 {
                    succ += remaining * up;
                    p_down += remaining * down;
                    remaining *= boom + stay;
                }
                let mut nonzero = 0;
                for (i, &n) in [du, ds, dd, db].iter().enumerate() {
                    if output[i][n as usize] == 0.0 {
                        uniq_subdists += 1;
                    }
                    if n != 0 {
                        nonzero += 1;
                    }
                    output[i][n as usize] += succ;
                }
                if nonzero > 1 && !joint.contains_key(&(du, ds, dd, db)) {
                    distr_adds += nonzero - 1;
                }
                // output[0][du as usize] += succ;
                // output[1][ds as usize] += succ;
                // output[2][dd as usize] += succ;
                // output[3][db as usize] += succ;
                merge_or_insert(&mut joint, (du, ds, dd, db), succ);
                total_prob += succ;
                update(states, (p_down, du, ds, dd, db, false));
            } else {
                update(states, (p * downup, du+1, ds, dd, db, true));
                update(states, (p * downstay, du, ds+1, dd, db, true));
                update(states, (p * downdown, du, ds, dd+1, db, true));
                update(states, (p * downboom, du, ds, dd, db+1, true));
                // dbg!((p, p*downdown, du, ds, dd+1, db));
                // assert!(p*downdown <= (down * downdown).powi((dd+1) as i32));
            }
        }
        let x = joint.get(&(0,0,0,0));
        dbg!(x);
        dbg!(output[0][0] * output[1][0] * output[2][0] * output[3][0]);
        dbg!(joint.len());
        dbg!(uniq_subdists);
        dbg!(distr_adds);
        output
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

    pub fn expected_cost(&self) -> u64 {
        let sum: f64 = self.dist.iter().map(|(c, p)| (*c as f64)*p).sum();
        let sum = (UNIT as f64) * sum;
        sum as u64
    }
}
