use std::hash::{BuildHasherDefault, Hash};

use indexmap::IndexMap;
use noisy_float::prelude::*;
use rustc_hash::FxHasher;

type F = R64;

pub struct Prio<K: Hash+Eq, V: Copy + Ord> {
    pub all_states: IndexMap<K, V, BuildHasherDefault<FxHasher>>,
    total_prob: F,
    pub prob: Option<fn(&mut V) -> &mut F>,
}

impl<K, V> Prio<K, V> where
    K: Hash + Eq, V: Copy + Ord
{
    pub fn new() -> Self {
        Self {
            all_states: Default::default(),
            total_prob: Default::default(),
            prob: Default::default()
        }
    }

    fn get(&self, idx: usize) -> V {
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

    fn heap_pop(&mut self) -> (K, V) {
        let popped = self.all_states.swap_remove_index(0);
        self.sift_down();
        unsafe {popped.unwrap_unchecked()}
    }

    pub fn push(&mut self, key: K, mut val: V) {
        let pf = &self.prob.as_ref().unwrap();
        let prob = pf(&mut val);
        self.total_prob += *prob;
        // merge two different probability paths of arriving at the same state
        let entry = self.all_states.entry(key);
        let idx = entry.index();
        entry.and_modify(|old| *pf(old) += *prob)
            .or_insert(val);
        self.sift_up(idx);
    }

    pub fn pop(&mut self) -> (K, V) {
        let (key, mut val) = self.heap_pop();
        self.total_prob -= *(self.prob.as_ref().unwrap())(&mut val);
        (key, val)
    }
}
