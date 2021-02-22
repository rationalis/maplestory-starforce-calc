use crate::binning::*;

use std::hash::Hash;
use std::ops::AddAssign;

use lazy_static::*;
use noisy_float::prelude::*;
use rustc_hash::FxHashMap;

pub type F = R64;

pub fn pp<T: Into<f64>>(num: T) -> String {
    use format_num::format_num;
    format_num!("0.3s", num)
}

pub fn f<T: Into<f64>>(f: T) -> F {
    r64(f.into())
}

pub fn g(f: F) -> f64 {
    f.raw()
}

pub type Meso = F;
pub type Star = u8;
pub type BinId = u16;

// pub const UNIT: Meso = 100_000;
pub const STAR_LIMIT: Star = 22;
pub const PROB_COUNT: usize = (STAR_LIMIT - 10) as usize;
pub const NUM_LEVELS: usize = 4;
pub const LEVELS: [i32; NUM_LEVELS] = [140, 150, 160, 200];

pub const MAX_BOOMS: usize = 15;
pub const MAX_DOWNS: usize = 64;

pub const DIST_THRESHOLD: f64 = 1e-6;
pub const IDENT_BINS: usize = 100;
pub const BASE_BIN: f64 = 1.62069;
pub const NUM_BINS: usize = 1300;
pub const BIN_EXP: f64 = 1.01;

pub const PROBS_F64: [[f64; 4]; PROB_COUNT] = [
    [0.5, 0.5, 0., 0.],
    [0.45, 0., 0.55, 0.],
    [0.4, 0., 0.594, 0.006],
    [0.35, 0., 0.637, 0.013],
    [0.3, 0., 0.686, 0.014],
    [0.3, 0.679, 0., 0.021], // 15
    [0.3, 0., 0.679, 0.021], // 16
    [0.3, 0., 0.679, 0.021], // 17
    [0.3, 0., 0.672, 0.028], // 18
    [0.3, 0., 0.672, 0.028], // 19
    [0.3, 0.63, 0., 0.07],   // 20
    [0.3, 0., 0.63, 0.07],   // 21
];

lazy_static! {
    pub static ref PROB_CUTOFF: F = f(1e-12);
    pub static ref ONE: F = f(1.0 - 0.5f64.powf(52.));
    pub static ref PROBS: [[F; 4]; PROB_COUNT] = {
        let mut probs: [[F; 4]; PROB_COUNT] = Default::default();
        for i in 0..12 {
            for j in 0..4 {
                probs[i][j] = f(PROBS_F64[i][j]);
            }
        }
        probs
    };
    pub static ref COST: [Meso; STAR_LIMIT as usize] = {
        let mut costs = [f(0); STAR_LIMIT as usize];
        for star in 10u8..22 {
            costs[star as usize] = cost_lv_indep(star);
        }
        costs
    };
    pub static ref BINS: [Meso; NUM_BINS] = {
        let mut bins = [f(0); NUM_BINS];
        for i in 0..NUM_BINS {
            if i <= IDENT_BINS {
                bins[i] = f(i as i32);
            } else {
                let frac: f64 = BIN_EXP.powi(i as i32 - IDENT_BINS as i32 - 1);
                bins[i] = f(BASE_BIN * frac);
            }
        }
        bins
    };
    pub static ref BINS_D: [Meso; NUM_BINS] = {
        let mut bins_d = [f(0); NUM_BINS];
        for i in 1..NUM_BINS {
            bins_d[i] = (BINS[i] + BINS[i - 1]) / 2.0;
        }
        bins_d
    };
    pub static ref BIN_SUMS: Vec<[BinId; NUM_BINS]> = {
        let mut bin_sums = Vec::with_capacity(NUM_BINS);
        for i in 0..NUM_BINS {
            let mut bin_sums_row = [0; NUM_BINS];
            for j in 0..NUM_BINS {
                bin_sums_row[j] = Bin::from_cost(BINS[i] + BINS[j]).raw();
            }
            bin_sums.push(bin_sums_row);
        }
        bin_sums
    };
}

// pub fn round(mesos: Meso, unit: i32) -> Meso {
//     (mesos + (unit / 2)) / unit * unit
// }

pub fn cost_lv_indep(star: Star) -> Meso {
    let star_base: Meso = f((star + 1) as f64);
    let star_pow: Meso = star_base.powf(f(2.7));
    let denom: Meso = f(match star {
        n if 10 <= n && n <= 14 => 400,
        n if 15 <= n && n <= 17 => 120,
        n if 18 <= n && n <= 19 => 110,
        n if 20 <= n && n <= 22 => 100,
        _ => panic!(),
    });
    let star_factor = star_pow / denom;
    println!("star factor at {} stars: {}", star, star_factor);
    star_factor
}

pub fn cost_conv_lv(c: Meso, level: i32) -> u64 {
    let level_factor = f(level.pow(3));
    let cost = level_factor * c;
    g(cost) as u64
    // let cost_approx = round(cost as i32, UNIT) / UNIT;
    // println!(
    //     "at level {} it costs {} mesos at {} stars, rounded to {}",
    //     level,
    //     pp(cost as u64),
    //     star,
    //     pp((round_bucket(cost_approx).1 * UNIT) as u64)
    // );
    // cost_approx
}

// pub fn cost(star: Star, level: i32) -> i32 {
//     let level_factor: f32 = level.pow(3) as f32;
//     let cost = 1000f32 + (level_factor as f32) * cost_lv_indep(star);
//     let cost_approx = round(cost as i32, UNIT) / UNIT;
//     println!(
//         "at level {} it costs {} mesos at {} stars, rounded to {}",
//         level,
//         pp(cost as u64),
//         star,
//         pp((round_bucket(cost_approx).1 * UNIT) as u64)
//     );
//     cost_approx
// }

pub fn merge_or_insert<K, V, D>(dist: &mut FxHashMap<K, V>, key: K, p: D)
where
    K: Eq + Hash,
    V: AddAssign<D>,
    D: Copy + Into<V>,
{
    dist.entry(key)
        .and_modify(|p0| *p0 += p)
        .or_insert(p.into());
}
