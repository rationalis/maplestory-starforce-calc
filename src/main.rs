#![feature(iter_map_while)]
#![feature(map_first_last)]
#![feature(option_result_unwrap_unchecked)]

use std::cmp::max;
use std::collections::{BTreeMap, BTreeSet};
use std::hash::{BuildHasherDefault, Hash, Hasher};
use std::time::SystemTime;

use indexmap::IndexMap;
use lazy_static::*;
use rustc_hash::{FxHasher, FxHashMap};
use noisy_float::prelude::*;
use serde::{Deserialize, Serialize};

use maplestory_calculator::prio::Prio;

type F = R64;
type Q = IndexMap<Key, F, BuildHasherDefault<FxHasher>>;
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

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Transition {
    Up = 0,
    Stay = 1,
    Down = 2,
    Boom = 3,
}

#[derive(Clone, Debug)]
struct Distr {
    dist: Vec<(Meso, f64)>
}

impl Distr {
    fn new(dist: Vec<(Meso, f64)>) -> Self {
        let mut map: FxHashMap<_, _> = Default::default();
        for (c, p) in dist.into_iter() {
            map.entry(round_bucket(c))
                .and_modify(|p0| *p0 += p)
                .or_insert(p);
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
        // // renormalize
        // self.dist.iter_mut().for_each(|(_, p)| *p /= total_prob);
    }

    fn constant(c: Meso) -> Self {
        Self {
            dist: vec![(c, 1.0)],
        }
    }

    fn zero() -> Self {
        Self::constant(0)
    }

    fn binom(n: i32, p:f64) -> Self {
        // println!("binom {} {}", n, p);

        // let mut dist = FxHashMap::default();
        let mut dist: [f64; 64] = [0.0; 64];
        let mut coeff: i64 = 1;
        let p2 = 1.0 - p;
        for k in 0..=n {
            if k > 0 {
                coeff *= (n - k + 1) as i64;
                coeff /= k as i64;
            }
            // if n >= 38 {
            //     println!("{} {} {}", coeff, p.powi(k), p2.powi(n-k));
            // }
            // dist.insert(k, (coeff as f64) * p.powi(k) * p2.powi(n-k));
            let p = coeff as f64 * p.powi(k) * p2.powi(n-k);
            dist[k as usize] = p;
        }
        // if n == 0 {
        //     let mut D: Vec<_> = dist.clone().into_iter().collect();
        //     D.sort_unstable_by_key(|(n,_)| *n);
        //     println!("{:?}", D);
        // }
        Self::new(dist.into_iter().enumerate().map(|(n, p)| (n as i32, *p)).collect())
    }

    fn geom(p: f64) -> Self {
        let mut dist = Vec::new();
        let mut S = 1.0;
        let mut i = 1;
        while S > 1e-6 {
            dist.push((i, S*p));
            i += 1;
            S -= S * p;
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
        // let update = |map: &mut Prio<State>, mut elem: State| {
        //     if elem.0 == 0. {
        //         return;
        //     }
        //     if map.contains(&elem) {
        //         map.remove(&elem);
        //         elem.0 *= 2.;
        //     }
        //     map.insert(elem);
        // };
        states.push((0,0,0,0,true), (f(1.0), 0));
        let mut output: [[f64; MAX_DOWN]; 4] = [[0.0; MAX_DOWN]; 4];
        let mut total_prob = 0.0;
        let mut i = 0;
        while total_prob < 1.0 - 1e-6 {
            i += 1;
            // if i < 100 {
            //     println!("{:?}", states)
            // } else { panic!() }
            if i % 100_000 == 0 {
                println!("{} {}", states.all_states.len(), total_prob);
            }
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
                if *[du, ds, dd, db].iter().max().unwrap() as usize == MAX_DOWN {
                    dbg!((du, ds, dd, db, at_start, p, total_prob));
                }
                output[0][du as usize] += succ;
                output[1][ds as usize] += succ;
                output[2][dd as usize] += succ;
                output[3][db as usize] += succ;
                total_prob += succ;
                // dbg!((du, ds, dd, db, at_start, p_down, succ));
                // panic!();
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
        Self::new(successes.into_iter().enumerate().map(|(b, p)| (b as i32, *p)).collect())
    }

    fn booms2(succ_rate: f64, boom_rate: f64) -> Self {
        const MAX_BOOMS: usize = 15;
        let mut dist: [f64; MAX_BOOMS] = [0.0; MAX_BOOMS];
        let mut fails = Distr::geom(succ_rate);
        fails.shift(-1);
        for (n, p) in fails.dist.iter() {
            let mut boom_subdist = Distr::binom(*n, boom_rate / (1.0 - succ_rate));
            for (n2, p2) in boom_subdist.dist.iter() {
                if *n2 < (MAX_BOOMS as i32) {
                    dist[*n2 as usize] += p*p2;
                }
                // dist.entry(*n2)
                //     .and_modify(|p0| *p0 += p*p2)
                //     .or_insert(p*p2);
            }
        }
        Self::new(dist.into_iter().enumerate().map(|(b, p)| (b as i32, *p)).collect())
    }

    fn prob_scale(&mut self, p: f64) -> &mut Self {
        self.dist.iter_mut().for_each(|(_, p0)| *p0 *= p);
        self
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
        println!("mult distrs of size {} and {}", self.dist.len(), other.dist.len());
        let mut dist = FxHashMap::default();
        for (i, p) in self.dist.iter() {
            for (c, p2) in other.dist.iter() {
                dist.entry(round_bucket(i*c))
                    .and_modify(|p0| *p0 += p*p2)
                    .or_insert(p*p2);
            }
        }
        let dist = dist.into_iter().collect();
        println!("finish mult distrs of size {} and {}", self.dist.len(), other.dist.len());
        Self::new(dist)
    }

    fn add(&self, other: &Self) -> Self {
        println!("add distrs of size {} and {}", self.dist.len(), other.dist.len());
        let mut dist = FxHashMap::default();
        for (c, p) in self.dist.iter() {
            for (c2, p2) in other.dist.iter() {
                dist.entry(round_bucket(c+c2))
                    .and_modify(|p0| *p0 += p*p2)
                    .or_insert(p*p2);
            }
        }
        let dist: Vec<_> = dist.into_iter().collect();
        println!("finish add distrs of size {} and {}", self.dist.len(), other.dist.len());
        println!("add result size: {}", dist.len());
        Self::new(dist)
    }

    fn expected_cost(&self) -> f64 {
        self.dist.iter().map(|(c, p)| (*c as f64)*p).sum()
    }
}

#[derive(Default, Serialize, Deserialize)]
struct Distribution {
    dist: Vec<(Meso, F)>,
    full_dist: Vec<(Meso, F)>
}

impl Distribution {
    fn new (full_dist: Vec<(Meso, F)>) -> Self {
        let mut small_dist = BTreeMap::new();

        for &(c, p) in full_dist.iter() {
            debug_assert!(c == round_bucket(c));
            small_dist.entry(c)
                .and_modify(|old_prob| *old_prob += p)
                .or_insert(p);
        }

        let mut dist: Vec<_> = small_dist.into_iter().collect();
        dist.sort_by_key(|(_, p)| *p);
        dist.reverse();

        Self {dist, full_dist}
    }
}

#[derive(Default)]
struct TransitionTable {
    dists: [[[Distribution; 2]; (STAR_LIMIT+2) as usize]; (STAR_LIMIT+2) as usize],
    highest: [[Star; 2]; (STAR_LIMIT+2) as usize],
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct State {
    prob: F,
    star: Star,
    spent: Meso,
    downed: bool,
}

#[derive(Eq, PartialEq)]
struct Key {
    star: Star,
    spent: Meso,
    downed: bool
}

impl Hash for Key {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let k = ((self.star as i32) << 24) | self.spent | ((self.downed as i32) << 31);
        k.hash(state);
    }
}

impl State {
    fn new(star: Star, downed: bool) -> Self {
        Self {
            prob: *ONE,
            star,
            spent: 0,
            downed,
        }
    }

    fn transition(&self, t: Transition) -> Option<Self> {
        use Transition::*;
        let p = PROBS[(self.star - 10) as usize][t as usize];
        if p == 0. {
            return None;
        }
        let mut next_downed = false;
        let mut spent = self.spent + COST[self.star as usize];
        let next_star = match t {
            Up => self.star + 1,
            Stay => self.star,
            Down if self.downed => {
                spent += COST[(self.star - 1) as usize];
                self.star
            },
            Down => {
                next_downed = true;
                self.star - 1
            },
            Boom => 12
        };
        let spent = round_bucket(spent);
        Some(Self {
            prob: self.prob * p,
            star: next_star,
            spent,
            downed: next_downed,
        })
    }

    fn join(self, table: &TransitionTable) -> Result<impl Iterator<Item=Self> + '_, Self> {
        let end = table.highest[self.star as usize][self.downed as usize];
        if end == 0 {
            return Err(self);
        }
        let dist = &table.dists[self.star as usize][end as usize][self.downed as usize].dist;
        Ok(dist.iter()
           .map_while(move |(c,p)| {
               let next_p = self.prob * p;
               if next_p > *PROB_CUTOFF {
                   Some(Self {
                       prob: next_p,
                       spent: round_bucket(self.spent + c),
                       star: end,
                       downed: false
                   })
               } else {
                   None
               }
           }))
    }

    fn add_transitions(&self, states: &mut States, table: &mut TransitionTable) {
        use Transition::*;
        for &t in &[Up, Stay, Down, Boom] {
            let next_state = self.transition(t);
            if let Some(next_state) = next_state {
                let next_states = next_state.join(table);
                match next_states {
                    Ok(next_states) => {
                        for next_state in next_states {
                            states.push(next_state);
                        }
                    },
                    Err(next_state) => states.push(next_state)
                }
            }
        }
    }

    fn key(&self) -> Key {
        Key {
            star: self.star,
            spent: self.spent,
            downed: self.downed,
        }
    }
}

#[derive(Default)]
struct States {
    all_states: Q,
    successes: FxHashMap<Meso, F>,
    total_prob: F,
    target: Star,

    success_prob: F,
    success_push_counter: i64,
    push_counter: i64,
    pop_counter: i64,
    merge_counter: i64,
    min_prob: F,
    max_len: i32,
}

impl States {
    fn get(&self, idx: usize) -> F {
        unsafe {*self.all_states.get_index(idx).unwrap_unchecked().1}
    }

    fn sift_up(&mut self, mut idx: usize) {
        while idx != 0 {
            let parent = (idx - 1) / 2;
            if self.get(parent) >= self.get(idx) { return; }
            self.all_states.swap_indices(idx, parent);
            idx = parent;
        }
    }

    fn sift_down(&mut self) {
        let len = self.all_states.len();
        let mut idx = 0;
        let mut left_child = 2*idx + 1;
        let mut right_child = 2*idx + 2;
        while left_child < len {
            let mut max_val = self.get(left_child);
            let mut max_child = left_child;
            if right_child < len {
                let right_val = self.get(right_child);
                if right_val > max_val {
                    max_val = right_val;
                    max_child = right_child;
                }
            }
            if max_val > self.get(idx) {
                self.all_states.swap_indices(idx, max_child);
                idx = max_child;
                left_child = 2*idx + 1;
                right_child = 2*idx + 2;
            } else {
                break;
            }
        }
    }

    fn heap_pop(&mut self) -> (Key, F) {
        let popped = self.all_states.swap_remove_index(0);
        self.pop_counter += 1;
        self.sift_down();
        unsafe {popped.unwrap_unchecked()}
    }

    fn push(&mut self, state: State) {
        let State {prob, star, spent, downed:_} = state;
        if star == self.target {
            self.success_prob += prob;
            self.success_push_counter += 1;
            self.successes.entry(spent)
                .and_modify(|v| {*v += prob;})
                .or_insert(prob);
            return;
        }
        self.total_prob += prob;
        self.min_prob = F::min(self.min_prob, prob);
        // merge two different probability paths of arriving at the same state
        let ct = &mut self.merge_counter;
        let entry = self.all_states.entry(state.key());
        let idx = entry.index();
        entry.and_modify(|old| {*old += prob; *ct += 1;})
            .or_insert(prob);
        self.sift_up(idx);

        self.push_counter += 1;
        self.max_len = max(self.max_len, self.all_states.len() as i32);
    }

    fn pop(&mut self) -> State {
        let (Key { star, spent, downed }, prob) = self.heap_pop();
        self.total_prob -= prob;
        State {prob, star, spent, downed}
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

fn calculate2(level: i32) {
    let mut ct = 0i64;
    let start_time = SystemTime::now();
    let mut table = TransitionTable::default();
    for target in 11..STAR_LIMIT+1 {
        for start in (10..target).rev() {
            for &downed in &[false, true] {
                println!("Starting {} downed: {} -> {}", start, downed, target);
                let states = calculate(State::new(start, downed), target, level, &mut table);
                let t = start_time.elapsed().unwrap();
                ct += states.push_counter;
                println!("{} pushes in {} secs", ct, t.as_secs_f32());
                println!("{} pushes/s", ct as f32 / t.as_secs_f32());
                let dist = Distribution::new(states.successes.into_iter().collect());
                table.dists[start as usize][target as usize][downed as usize] = dist;
                table.highest[start as usize][downed as usize] = target;
            }
        }
    }
    let s = serde_json::to_string(&table.dists[10][22][0]).expect("failed tostr");
    std::fs::write("output.json", s).expect("failed file write");
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
            let mut base = Distr::geom(up);
            // cost of attempts at start->target
            let mut base2 = base.clone();
            base2.scale(COST[start as usize]);
            dist = dist.add(&base2);

            if down > 0. {
                let cost_below = COST[(start - 1) as usize];
                let [downup, _downstay, downdown, downboom] = PROBS_F64[(start - 11) as usize];
                // if we can go down then we need the number of failures
                if downdown == 0. {
                    let mut fails = base.clone();
                    fails.shift(-1);
                    let DS_cost = table.get(&(start-1, start)).unwrap();
                    dist = dist.add(&fails.product(&DS_cost));
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

fn calculate(start: State, target: Star, level: i32, table: &mut TransitionTable) -> States {
    let mut states = States {
        target,
        min_prob: *ONE,
        ..Default::default()};
    let mut last_total_prob = *ONE;
    let threshold = 1e-6;
    states.push(start);
    while states.total_prob > threshold {
        states.pop().add_transitions(&mut states, table);
        // if true {
        if cfg!(debug_assertions) {
            let mut progressed = false;
            for magnitude in 1.. {
                let oom = 0.1f64.powf(magnitude as f64);
                if last_total_prob > oom {
                    progressed = states.total_prob <= last_total_prob - f(oom);
                    break;
                }
                if oom < threshold {
                    break;
                }
            }
            if progressed {
                println!("{}", states.total_prob);
                last_total_prob = states.total_prob;
            }
        }
    }
    if true {
    // if cfg!(debug_assertions) {
        let mut expected_cost = 0.;
        for (c, p) in states.successes.iter() {
            let p: f64 = g(*p);
            expected_cost += p*(*c as f64);
        }
        println!("Expected cost: {}", (expected_cost * UNIT as f64) as i64);
        println!("Max cost: {}", (*states.successes.keys().max().unwrap() as f64 * UNIT as f64) as i64);
        println!("States popped: {}", states.pop_counter);
        println!("States pushed: {}, Merges: {}", states.push_counter, states.merge_counter);
        println!("Unmerged states pushed: {}", states.push_counter - states.merge_counter);
        println!("Success states pushed: {}", states.success_push_counter);
        println!("Vanished prob: {}", 1.0 - g(states.total_prob) - g(states.success_prob));
        println!("Smallest non-success prob path: {}", states.min_prob);
        println!("Max len: {}", states.max_len);
        println!("End len: {}", states.all_states.len());
        println!("Number of success paths: {}", states.successes.len());
    }
    states
}

fn main() {
    // calculate(10, 17, 160, &mut TransitionTable::default());
    // calculate2(160);
    calculate3(160);
}
