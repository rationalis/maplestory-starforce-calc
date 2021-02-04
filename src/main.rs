#![feature(map_first_last)]

use std::collections::{BTreeSet, HashMap};

use lazy_static::*;
use noisy_float::prelude::*;

const UNIT: i32 = 100_000;

const PROBS: [[f32; 4]; 7] = [
    [ 0.5, 0.5, 0., 0. ],
    [ 0.45, 0., 0.55, 0. ],
    [ 0.4, 0., 0.594, 0.006 ],
    [ 0.35, 0., 0.637, 0.013 ],
    [ 0.3, 0., 0.686, 0.014 ],
    [ 0.3, 0.679, 0., 0.021 ],
    [ 0.3, 0., 0.679, 0.021 ]
];

lazy_static! {
    static ref LEVEL: i32 = 160;
    static ref COST: [i32; 17] = {
        let mut cost = [0; 17];

        for i in 10..17 {
            let i = i as usize;
            cost[i] = crate::cost(i as i32, *LEVEL);
        }

        cost
    };
}

fn round(mesos: i32, unit: i32) -> i32 {
    (mesos + (unit / 2)) / unit * unit
}


#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Transition {
    Up = 0,
    Stay = 1,
    Down = 2,
    Boom = 3,
}

#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
struct State {
    prob: N32,
    star: i32,
    spent: i32,
    downed: bool,
}

impl State {
    fn new(star: i32) -> Self {
        Self {
            prob: n32(1.),
            star,
            spent: 0,
            downed: false,
        }
    }

    fn transition(&self, t: Transition) -> Option<Self> {
        use Transition::*;
        let p = PROBS[(self.star - 10) as usize][t as usize];
        if p == 0. {
            return None;
        }
        let mut spent = self.spent + COST[self.star as usize];
        let next_star = match t {
            Up => self.star + 1,
            Stay => self.star,
            Down if self.downed => self.star,
            Down => self.star - 1,
            Boom => 12
        };
        if self.downed && t == Transition::Down {
            spent += COST[(self.star - 1) as usize];
        }
        if spent > 50000 {
            spent = round(spent, 100);
        } else if spent > 5000 {
            spent = round(spent, 10);
        }
        Some(Self {
            prob: self.prob * p,
            star: next_star,
            spent,
            downed: t == Transition::Down && !self.downed,
        })
    }

    fn add_transitions(&self, states: &mut States) {
        use Transition::*;
        for t in &[Up, Stay, Down, Boom] {
            let next_state = self.transition(*t);
            if let Some(next_state) = next_state {
                states.push(next_state);
            }
        }
    }

    fn key(&self) -> (i32, i32, bool) {
        (self.star, self.spent, self.downed)
    }
}

#[derive(Default)]
struct States {
    state_to_prob: HashMap<(i32, i32, bool), N32>,
    all_states: BTreeSet<State>,
    successes: Vec<(N32, i32)>,
    total_prob: N32,
    target: i32,
    push_counter: i32,
    merge_counter: i32,
}

impl States {
    fn push(&mut self, mut state: State) {
        let State {prob, star, spent, downed:_} = state;
        if star == self.target {
            self.successes.push((prob, spent));
            return;
        } else {
            self.total_prob += prob;
        }
        // merge two different probability paths of arriving at the same state
        if let Some(&old_prob) = self.state_to_prob.get(&state.key()) {
            let new_prob = old_prob + prob;
            self.state_to_prob.insert(state.key(), new_prob);
            let mut copy: State = state;
            copy.prob = old_prob;
            state.prob = new_prob;
            self.all_states.remove(&copy);
            self.all_states.insert(state);
            self.merge_counter += 1;
        } else {
            self.state_to_prob.insert(state.key(), prob);
            self.all_states.insert(state);
        }
        self.push_counter += 1;
    }

    fn pop(&mut self) -> State {
        let state = self.all_states.pop_last().unwrap();
        self.state_to_prob.remove(&state.key());
        self.total_prob -= state.prob;
        state
    }
}

fn cost(star: i32, level: i32) -> i32 {
    let level_factor: f32 = level.pow(3) as f32;
    let star_factor: f32 = ((star+1) as f32).powf(2.7);
    let denom: f32 = match star {
        n if 10 <= n && n <= 14 => 400,
        n if 15 <= n && n <= 17 => 120,
        _ => panic!(),
    } as f32;
    let cost = 1000f32 + (level_factor as f32) * star_factor / denom;
    let cost_approx = round(cost as i32, UNIT) / UNIT;
    println!("costs {} mesos at {} stars, rounded to {}", cost, star, cost_approx * UNIT);
    cost_approx
}

fn calculate() {
    let mut states = States {target: 17, ..Default::default()};
    let mut last_total_prob = n32(1.0);
    let threshold = 1e-6;
    let init_state = State::new(10);
    states.push(init_state);
    while states.total_prob > threshold {
        states.pop().add_transitions(&mut states);
        let mut progressed = false;
        for magnitude in 2.. {
            let oom = 0.1f32.powf(magnitude as f32);
            if last_total_prob > oom {
                progressed = states.total_prob <= last_total_prob - oom;
                break;
            }
            if oom < threshold {
                break;
            }
        }
        if last_total_prob > 0.01 {
            states.total_prob <= last_total_prob - 0.01
        } else {
            states.total_prob <= last_total_prob - 0.001
        };
        if progressed {
            println!("{}", states.total_prob);
            last_total_prob = states.total_prob;
        }
    }
    let mut expected_cost = 0.;
    for (p, c) in states.successes {
        let p: f32 = p.into();
        expected_cost += p*(c as f32);
    }
    println!("Expected cost: {}", expected_cost * UNIT as f32);
    println!("States pushed: {}, Merges: {}", states.push_counter, states.merge_counter);
}

fn main() {
    calculate();
}
