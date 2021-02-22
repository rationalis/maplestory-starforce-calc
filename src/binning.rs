use crate::consts::*;

use rustc_hash::FxHashMap;
use nanorand::{RNG, WyRand};

pub fn round_bucket_impl(bins: &[Meso], mesos: Meso) -> (BinId, Meso) {
    let c = mesos;
    let res = bins.binary_search(&c);
    let idx = match res {
        Ok(found) => found,
        Err(insertion) => {
            let mut closest = None;
            let mut closest_diff = None;
            if insertion < bins.len() {
                closest = Some(insertion);
                closest_diff = Some(g(bins[insertion] - c).abs());
            }
            if insertion > 0 {
                if closest.is_none() || g(bins[insertion - 1] - c).abs() < closest_diff.unwrap() {
                    closest = Some(insertion - 1);
                }
            }
            closest.unwrap()
        }
    };
    (idx as BinId, bins[idx])
}

#[repr(transparent)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd, Hash)]
pub struct Bin {
    bin_id: BinId,
}

impl Bin {
    #[inline(always)]
    pub fn from_id(bin_id: BinId) -> Self {
        Self { bin_id }
    }

    #[inline(always)]
    pub fn from_int(n: i32) -> Self {
        Self { bin_id: n as u16 }
    }

    #[inline(always)]
    pub fn from_cost(c: Meso) -> Self {
        let idx: BinId = round_bucket_impl(&BINS[IDENT_BINS + 1..], c).0;
        Self { bin_id: idx + IDENT_BINS as BinId + 1 }
    }

    #[inline(always)]
    pub fn raw(&self) -> BinId {
        self.bin_id
    }

    #[inline(always)]
    pub fn unbin(&self) -> Meso {
        BINS[self.bin_id as usize]
    }
}

pub fn trajectory(start: Star, target: Star, rng: &mut WyRand,
                  counters: &mut FxHashMap<Meso, i32>) {
    let mut cost = f(0.0);
    let mut star = start;
    while star != target {
        // add cost, update counters
        cost += COST[star as usize];
        merge_or_insert(counters, cost, 1);

        // simulate one transition
        let [up, _, down, boom] = PROBS_F64[(star - 10) as usize];
        let rand: f64 = rng.generate_range::<u32>(0, 1_000_000) as f64 * 1e-6;
        if rand < up {
            star += 1;
        } else if rand < up + down {
            star -= 1;
        } else if rand < up + down + boom {
            star = 12;
        }
    }
}

pub fn bins() -> Vec<Meso> {
    const THRESH: f64 = 0.005;
    let mut counters = Default::default();
    let mut rng = WyRand::new_seed(0);
    for start in 10..=21 {
        for target in start+1..=22 {
            trajectory(start, target, &mut rng, &mut counters);
        }
    }
    for _ in 0..100 {
        trajectory(10, 22, &mut rng, &mut counters);
    }
    dbg!(counters.len());
    let mut list: Vec<_> = counters.into_iter().collect();
    list.sort_unstable_by_key(|(cost, _)| *cost);
    let mut means = Vec::new();
    let mut i = 0;
    while i < list.len() {
        let lower = list[i].0;
        let mut cs = list[i].0 * f(list[i].1); // current sum
        let mut cc = list[i].1; // current count
        i += 1;
        while i < list.len() {
            let (next_cost, next_count) = list[i];
            let mult_cost = next_cost * f(next_count);
            let next_mean = (cs + mult_cost) / f(cc + next_count);
            let small_enough = next_mean <= lower * (1.0 + THRESH);
            let big_enough = next_mean >= next_cost * (1.0 - THRESH);
            if small_enough && big_enough {
                cs += mult_cost;
                cc += next_count;
                i += 1;
            } else {
                break;
            }
        }
        let mean = cs / f(cc);
        means.push(mean);
    }
    dbg!(means.len());
    means
}

