#![feature(map_first_last)]
#![feature(option_result_unwrap_unchecked)]

use std::collections::{BTreeMap, HashMap};

use fixed::types::*;
use lazy_static::*;
use indexmap::IndexMap;

type F = U1F63;
type Q = IndexMap<Key, F, std::hash::BuildHasherDefault<rustc_hash::FxHasher>>;
type Meso = i32;
type Star = u8;

const UNIT: Meso = 100_000;
const STAR_LIMIT: Star = 22;
const PROB_COUNT: usize = (STAR_LIMIT - 10) as usize;

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
    // TODO: add level-dependent cost tables
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
    if c > 10000 {
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
            let c = round_bucket(c);
            if let Some(&old_prob) = small_dist.get(&c) {
                small_dist.insert(c, old_prob + p);
            } else {
                small_dist.insert(c, p);
            }
        }

        let dist = small_dist.into_iter().collect();

        Self {dist, full_dist}
    }
}

#[derive(Default)]
struct TransitionTable {
    dists: HashMap<(Star, bool), BTreeMap<Star, Distribution>>
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
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

impl std::hash::Hash for Key {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let k = ((self.star as i32) << 24) | self.spent | ((self.downed as i32) << 31);
        k.hash(state);
    }
}

impl State {
    fn new(star: Star, downed: bool) -> Self {
        Self {
            prob: F::from_num(1.),
            star,
            spent: 0,
            downed,
        }
    }

    #[inline(always)]
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

    fn join(self, table: &TransitionTable) -> Box<dyn Iterator<Item=State> + '_> {
        let mut dist = None;
        if let Some(starting) = table.dists.get(&(self.star, self.downed)) {
            if let Some(ending) = starting.last_key_value() {
                dist = Some(ending);
            }
        }
        if dist.is_none() {
            return Box::new(std::iter::once(self));
        }
        let (end, dist) = dist.unwrap();
        let end = *end;
        let dist = &dist.dist;
        Box::new(dist.iter().map(move |(c,p)| Self {
            prob: self.prob * p,
            spent: self.spent + c,
            star: end,
            downed: false,
        }).filter(|s| s.prob > 0))
    }

    fn add_transitions(&self, states: &mut States, table: &mut TransitionTable) {
        use Transition::*;
        for &t in &[Up, Stay, Down, Boom] {
            let next_state = self.transition(t);
            if let Some(next_state) = next_state {
                let next_states = next_state.join(table);
                for next_state in next_states {
                    states.push(next_state);
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

    push_counter: i64,
    merge_counter: i64,
    min_prob: F,
    max_len: i32,
}

impl States {
    fn get(&self, idx: usize) -> F {
        unsafe {*self.all_states.get_index(idx).unwrap_unchecked().1}
    }

    fn sift_up(&mut self, mut idx: usize) {
        let mut parent = (idx - 1) / 2;
        while idx != 0 && self.get(parent) < self.get(idx) {
            self.all_states.swap_indices(idx, parent);
            idx = parent;
            parent = (idx - 1) / 2;
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
        self.sift_down();
        unsafe {popped.unwrap_unchecked()}
    }

    fn push(&mut self, state: State) {
        let State {prob, star, spent, downed:_} = state;
        if star == self.target {
            if let Some(old_prob) = self.successes.get_mut(&spent) {
                *old_prob += prob;
            } else {
                self.successes.insert(spent, prob);
            }
            return;
        } else {
            self.total_prob += prob;
            self.min_prob = F::min(self.min_prob, prob);
        }
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
    let mut table = TransitionTable::default();
    for target in 11..STAR_LIMIT {
        for start in (10..target).rev() {
            for &downed in &[false, true] {
                println!("Starting {} downed: {} -> {}", start, downed, target);
                let states = calculate(State::new(start, downed), target, level, &mut table);
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
        min_prob: F::from_num(1.),
        ..Default::default()};
    let mut last_total_prob = F::from_num(1.0);
    let threshold = if target <= 17 { 1e-5 } else { 1e-4 };
    states.push(start);
    while states.total_prob > threshold {
        states.pop().add_transitions(&mut states, table);
        if true {
        // if cfg!(debug_assertions) {
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
    let mut expected_cost = 0.;
    for (c, p) in states.successes.iter() {
        let p: f64 = p.to_num(); //p.into();
        expected_cost += p*(*c as f64);
    }
    println!("Expected cost: {}", (expected_cost * UNIT as f64) as i32);
    println!("States pushed: {}, Merges: {}", states.push_counter, states.merge_counter);
    println!("Unmerged states pushed: {}", states.push_counter - states.merge_counter);
    println!("Smallest prob path: {}", states.min_prob);
    println!("Max len: {}", states.max_len);
    println!("End len: {}", states.all_states.len());
    println!("Number of success paths: {}", states.successes.len());
    states
}

fn main() {
    // calculate(10, 17, 160, &mut TransitionTable::default());
    calculate2(160);
}
