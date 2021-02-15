#![feature(iter_map_while)]
#![feature(map_first_last)]
#![feature(option_result_unwrap_unchecked)]

use std::time::SystemTime;

use lazy_static::*;
use rustc_hash::FxHashMap;
use noisy_float::prelude::*;
use serde::{Deserialize, Serialize};

use maplestory_calculator::prio::Prio;

type F = R64;
type Meso = i32;
type Star = u8;

const UNIT: Meso = 100_000;
const STAR_LIMIT: Star = 22;
const PROB_COUNT: usize = (STAR_LIMIT - 10) as usize;

const MAX_DOWN: usize = 64;

// TODO: also calculate # booms (independently)
// TODO: add level-dependent cost tables
// TODO: add options for star-catching, event discounts
// TODO: command-line options
// TODO: save outputs to a file
// TODO: plot outputs
// TODO: write README

const PROBS_F64: [[f64; 4]; PROB_COUNT] = [
    [ 0.5, 0.5, 0., 0. ],
    [ 0.45, 0., 0.55, 0. ],
    [ 0.4, 0., 0.594, 0.006 ],
    [ 0.35, 0., 0.637, 0.013 ],
    [ 0.3, 0., 0.686, 0.014 ],
    [ 0.3, 0.679, 0., 0.021 ], // 15
    [ 0.3, 0., 0.679, 0.021 ], // 16
    [ 0.3, 0., 0.679, 0.021 ], // 17
    [ 0.3, 0., 0.672, 0.028 ], // 18
    [ 0.3, 0., 0.672, 0.028 ], // 19
    [ 0.3, 0.63, 0., 0.07 ], // 20
    [ 0.3, 0., 0.63, 0.07 ], // 21
];

lazy_static! {
    static ref PROB_CUTOFF: F = f(1e-12);
    static ref ONE: F = f(1.0 - 0.5f64.powf(52.));
    static ref PROBS: [[F; 4]; PROB_COUNT] = {
        let mut probs: [[F; 4]; PROB_COUNT] = Default::default();
        for i in 0..12 {
            for j in 0..4 {
                probs[i][j] = f(PROBS_F64[i][j]);
            }
        }
        probs
    };
    static ref LEVEL: i32 = 160;
    static ref COST: [Meso; STAR_LIMIT as usize] = {
        let mut cost = [0; STAR_LIMIT as usize];

        for i in 10..22 {
            let i = i as usize;
            cost[i] = crate::cost(i as u8, *LEVEL);
        }

        cost
    };
}

fn f(f: f64) -> F {
    r64(f)
}

fn g(f: F) -> f64 {
    f.raw()
}

fn slice_to_distr(slice: &[f64]) -> Distr {
    Distr::new(slice.into_iter().enumerate().map(|(n, p)| (n as i32, *p)).collect())
}

fn merge_or_insert(dist: &mut FxHashMap<Meso, f64>, key: Meso, p: f64) {
    dist.entry(key)
        .and_modify(|p0| *p0 += p)
        .or_insert(p);
}

fn round(mesos: Meso, unit: i32) -> Meso {
    (mesos + (unit / 2)) / unit * unit
}

fn round_bucket(mesos: Meso) -> Meso {
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

#[derive(Clone, Debug)]
struct Distr {
    dist: Vec<(Meso, f64)>
}

impl Distr {
    fn new(dist: Vec<(Meso, f64)>) -> Self {
        let mut map: FxHashMap<_, _> = Default::default();
        for (c, p) in dist.into_iter() {
            merge_or_insert(&mut map, round_bucket(c), p);
        }
        let dist: Vec<_> = map.into_iter().collect();
        let mut res = Self { dist };
        res.truncate(1e-6);
        res
    }

    fn truncate(&mut self, threshold: f64) {
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
        self.dist.push((last_key, 1.0 - total_prob));
    }

    fn constant(c: Meso) -> Self {
        Self {
            dist: vec![(c, 1.0)],
        }
    }

    fn zero() -> Self {
        Self::constant(0)
    }

    fn geom(p: f64) -> Self {
        let mut dist = Vec::new();
        let mut remaining = 1.0;
        let mut i = 1;
        while remaining > 1e-6 {
            dist.push((i, remaining*p));
            i += 1;
            remaining -= remaining * p;
        }
        Self {
            dist
        }
    }

    fn downs(start: Star) -> [[f64; MAX_DOWN]; 4] {
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
            map.push((du, ds, dd, db, at_start), (p, -((du + ds + dd + db) as i32)));
        };
        states.push((0,0,0,0,true), (f(1.0), 0));
        let mut output: [[f64; MAX_DOWN]; 4] = [[0.0; MAX_DOWN]; 4];
        let mut total_prob = 0.0;
        while total_prob < 1.0 - 1e-6 {
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
                output[0][du as usize] += succ;
                output[1][ds as usize] += succ;
                output[2][dd as usize] += succ;
                output[3][db as usize] += succ;
                total_prob += succ;
                update(states, (p_down, du, ds, dd, db, false));
            } else {
                update(states, (p * downup, du+1, ds, dd, db, true));
                update(states, (p * downstay, du, ds+1, dd, db, true));
                update(states, (p * downdown, du, ds, dd+1, db, true));
                update(states, (p * downboom, du, ds, dd, db+1, true));
            }
        }
        output
    }

    /// This method is actually faster by something like 25x compared to summing
    /// the probabilities over the compound binomial distribution.
    fn booms(succ_rate: f64, boom_rate: f64) -> Self {
        const MAX_BOOMS: usize = 15;
        let mut states: [f64; MAX_BOOMS] = [0.0; MAX_BOOMS];
        let mut successes: [f64; MAX_BOOMS] = [0.0; MAX_BOOMS];
        states[0] = 1.0;
        let mut total_prob = 0.0;
        while total_prob < 1.0 - 1e-9 {
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

    fn shift(&mut self, c: Meso) -> &mut Self {
        self.dist.iter_mut().for_each(|(c0, _)| *c0 += c);
        self
    }

    fn scale(&mut self, n: Meso) -> &mut Self {
        self.dist.iter_mut().for_each(|(c0, _)| *c0 *= n);
        self
    }

    fn product(&self, other: &Self) -> Self {
        let mut dist = FxHashMap::default();
        for (i, p) in self.dist.iter() {
            for (c, p2) in other.dist.iter() {
                merge_or_insert(&mut dist, round_bucket(i*c), p*p2);
            }
        }
        let dist = dist.into_iter().collect();
        Self::new(dist)
    }

    fn add(&self, other: &Self) -> Self {
        let mut dist = FxHashMap::default();
        for (c, p) in self.dist.iter() {
            for (c2, p2) in other.dist.iter() {
                merge_or_insert(&mut dist, round_bucket(c+c2), p*p2);
            }
        }
        let dist: Vec<_> = dist.into_iter().collect();
        Self::new(dist)
    }

    fn expected_cost(&self) -> f64 {
        self.dist.iter().map(|(c, p)| (*c as f64)*p).sum()
    }
}

fn cost(star: Star, level: i32) -> i32 {
    let level_factor: f32 = level.pow(3) as f32;
    let star_factor: f32 = ((star+1) as f32).powf(2.7);
    let denom: f32 = match star {
        n if 10 <= n && n <= 14 => 400,
        n if 15 <= n && n <= 17 => 120,
        n if 18 <= n && n <= 19 => 110,
        n if 20 <= n && n <= 22 => 100,
        _ => panic!(),
    } as f32;
    let cost = 1000f32 + (level_factor as f32) * star_factor / denom;
    let cost_approx = round(cost as i32, UNIT) / UNIT;
    println!("costs {} mesos at {} stars, rounded to {}", cost, star, cost_approx * UNIT);
    cost_approx
}

fn calculate3(level: i32) {
    let mut table: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    table.insert((12, 12), Distr::zero());
    for target in 11..STAR_LIMIT+1 {
        for start in (10..target).rev() {
            println!("{} -> {}", start, target);
            if start != target - 1 {
                let dist1 = table.get(&(start, start+1)).unwrap();
                let dist2 = table.get(&(start+1, target)).unwrap();
                let dist = dist1.add(dist2);
                println!("storing dist of size: {}", dist.dist.len());
                println!("Expected cost: {}", (dist.expected_cost() * UNIT as f64) as i64);
                table.insert((start, target), dist);
                continue;
            }
            let [up, _stay, down, boom] = PROBS_F64[(start - 10) as usize];
            if down == 0. && boom == 0. {
                let mut dist = Distr::geom(up);
                dist.scale(COST[start as usize]);
                println!("storing dist of size: {}", dist.dist.len());
                println!("Expected cost: {}", (dist.expected_cost() * UNIT as f64) as i64);
                table.insert((start, target), dist);
                continue;
            }
            let mut dist = Distr::zero();
            let base = Distr::geom(up);
            // cost of attempts at start->target
            let mut base2 = base.clone();
            base2.scale(COST[start as usize]);
            dist = dist.add(&base2);

            if down > 0. {
                let cost_below = COST[(start - 1) as usize];
                let [_, _, downdown, downboom] = PROBS_F64[(start - 11) as usize];
                // if we can go down then we need the number of failures
                if downdown == 0. {
                    let mut fails = base.clone();
                    fails.shift(-1);
                    let ds_cost = table.get(&(start-1, start)).unwrap();
                    dist = dist.add(&fails.product(&ds_cost));
                } else {
                    let cost_two_below = COST[(start - 2) as usize];
                    // println!("Starting chance time sim");
                    let [downups, downstays, downdowns, downbooms] = Distr::downs(start);
                    let downups = slice_to_distr(&downups);
                    let downstays = slice_to_distr(&downstays);
                    let downdowns = slice_to_distr(&downdowns);
                    let downbooms = slice_to_distr(&downbooms);
                    // dbg!(downdowns);
                    // println!("Finished chance time sim");
                    let mut du_cost = downups;
                    du_cost.scale(cost_below);
                    dist = dist.add(&du_cost);
                    let mut ds_cost = table.get(&(start - 1, start)).unwrap().clone();
                    ds_cost.shift(cost_below);
                    ds_cost = ds_cost.product(&downstays);
                    dist = dist.add(&ds_cost);
                    let mut dd_cost = table.get(&(start - 1, start)).unwrap().clone();
                    dd_cost.shift(cost_below).shift(cost_two_below);
                    dd_cost = dd_cost.product(&downdowns);
                    dist = dist.add(&dd_cost);
                    if downboom > 0. {
                        let mut db_cost = table.get(&(12, start)).unwrap().clone();
                        db_cost.shift(cost_below);
                        db_cost = db_cost.product(&downbooms);
                        dist = dist.add(&db_cost);
                    }
                }
            }
            if boom > 0. {
                let booms = Distr::booms(up, boom);
                let boom_cost = table.get(&(12, start)).unwrap();
                dist = dist.add(&booms.product(boom_cost));
            }
            println!("storing dist of size: {}", dist.dist.len());
            println!("Expected cost: {}", (dist.expected_cost() * UNIT as f64) as i64);
            table.insert((start, target), dist);
        }
    }
}

fn main() {
    calculate3(160);
}
