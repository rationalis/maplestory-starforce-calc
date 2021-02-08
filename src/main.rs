#![feature(iter_map_while)]
#![feature(map_first_last)]
#![feature(option_result_unwrap_unchecked)]

use std::collections::{BTreeMap, HashMap};
use std::hash::{BuildHasherDefault, Hash, Hasher};
use std::time::SystemTime;

use fixed::types::U0F64;
use indexmap::IndexMap;
use lazy_static::*;
use rustc_hash::FxHasher;

type F = U0F64;
type Q = IndexMap<Key, F, BuildHasherDefault<FxHasher>>;
type Meso = i32;
type Star = u8;

const UNIT: Meso = 100_000;
const STAR_LIMIT: Star = 22;
const PROB_COUNT: usize = (STAR_LIMIT - 10) as usize;

// TODO: also calculate # booms (independently)
// TODO: add level-dependent cost tables
// TODO: add options for star-catching, event discounts
// TODO: command-line options
// TODO: save outputs to a file
// TODO: plot outputs
// TODO: write README

const PROBS_F32: [[f32; 4]; PROB_COUNT] = [
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
    static ref PROB_CUTOFF: F = F::from_num(1e-14);
    static ref ONE: F = F::from_num(1.0 - 0.5f64.powf(52.));
    static ref PROBS: [[F; 4]; PROB_COUNT] = {
        let mut probs: [[F; 4]; PROB_COUNT] = Default::default();
        for i in 0..12 {
            for j in 0..4 {
                probs[i][j] = F::from_num(PROBS_F32[i][j]);
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

#[derive(Default)]
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
    dists: HashMap<(Star, bool), BTreeMap<Star, Distribution>>
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
        let mut dist = None;
        if let Some(starting) = table.dists.get(&(self.star, self.downed)) {
            if let Some(ending) = starting.last_key_value() {
                dist = Some(ending);
            }
        }
        if dist.is_none() {
            return Err(self);
        }
        let (end, dist) = dist.unwrap();
        let end = *end;
        let dist = &dist.dist;
        Ok(dist.iter()
           .map_while(move |(c,p)| {
               let next_p = self.prob * p;
               if next_p.to_bits() > PROB_CUTOFF.to_bits() {
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
    successes: HashMap<Meso, F>,
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
        self.max_len = std::cmp::max(self.max_len, self.all_states.len() as i32);
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
                if let Some(starting) = table.dists.get_mut(&(start, downed)) {
                    starting.insert(target, dist);
                } else {
                    let mut map = BTreeMap::new();
                    map.insert(target, dist);
                    table.dists.insert((start, downed), map);
                }
            }
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
                let oom = 0.1f32.powf(magnitude as f32);
                if last_total_prob > oom {
                    progressed = states.total_prob <= last_total_prob - F::from_num(oom);
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
    if cfg!(debug_assertions) {
        let mut expected_cost = 0.;
        for (c, p) in states.successes.iter() {
            let p: f64 = p.to_num(); //p.into();
            expected_cost += p*(*c as f64);
        }
        println!("Expected cost: {}", (expected_cost * UNIT as f64) as i64);
        println!("Max cost: {}", (*states.successes.keys().max().unwrap() as f64 * UNIT as f64) as i64);
        println!("States popped: {}", states.pop_counter);
        println!("States pushed: {}, Merges: {}", states.push_counter, states.merge_counter);
        println!("Unmerged states pushed: {}", states.push_counter - states.merge_counter);
        println!("Success states pushed: {}", states.success_push_counter);
        println!("Vanished prob: {}", 1.0 - states.total_prob.to_num::<f64>() - states.success_prob.to_num::<f64>());
        println!("Smallest non-success prob path: {}", states.min_prob);
        println!("Max len: {}", states.max_len);
        println!("End len: {}", states.all_states.len());
        println!("Number of success paths: {}", states.successes.len());
    }
    states
}

fn main() {
    // calculate(10, 17, 160, &mut TransitionTable::default());
    calculate2(160);
}
