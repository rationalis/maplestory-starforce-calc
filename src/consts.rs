use crate::distr::round_bucket;

use lazy_static::*;
use noisy_float::prelude::*;

pub type F = R64;

pub fn pp(num: u64) -> String {
    use format_num::format_num;
    format_num!("0.3s", num as f64)
}

pub fn f(f: f64) -> F {
    r64(f)
}

pub fn g(f: F) -> f64 {
    f.raw()
}

pub type Meso = f64;
pub type Star = u8;

pub const BASE_STAR_FACTOR: f32 = 11f32.powf(2.7) / 400f32;

pub const STAR_LIMIT: Star = 22;
pub const PROB_COUNT: usize = (STAR_LIMIT - 10) as usize;
pub const NUM_LEVELS: usize = 4;
pub const LEVELS: [i32; NUM_LEVELS] = [140, 150, 160, 200];

pub const MAX_BOOMS: usize = 15;
pub const MAX_DOWNS: usize = 64;

pub const DIST_THRESHOLD: f64 = 1e-6;
pub const NUM_BINS: usize = 2048;
pub const BIN_EXP: f64 = 1.0085;

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
    pub static ref COST: [[Meso; STAR_LIMIT as usize]; NUM_LEVELS] = {
        let mut costs = [[0f32.into(); STAR_LIMIT as usize]; NUM_LEVELS];

        for (i, &level) in LEVELS.iter().enumerate() {
            for star in 10u8..22 {
                //costs[i as usize][star as usize] = cost(star, level);
                costs[i as usize][star as usize] = cost(star, level);
            }
        }

        costs
    };
    pub static ref BINS: [Meso; NUM_BINS] = {
        let mut bins = [0f32.into(); NUM_BINS];
        for i in 0..NUM_BINS {
            let frac: f64 = BIN_EXP.powi(i as i32);
            bins[i] = ((frac as f32) * BASE_STAR_FACTOR).into();
        }
        bins
    };
    pub static ref BIN_SUMS: Vec<[u16; NUM_BINS]> = {
        let mut bin_sums = Vec::with_capacity(NUM_BINS);
        for i in 0..NUM_BINS {
            let mut bin_sums_row = [0; NUM_BINS];
            for j in 0..NUM_BINS {
                bin_sums_row[j] = round_bucket(BINS[i] + BINS[j]).0 as u16;
            }
            bin_sums.push(bin_sums_row);
        }
        bin_sums
    };
}

// pub fn round(mesos: Meso, unit: i32) -> Meso {
//     (mesos + (unit / 2)) / unit * unit
// }

pub fn cost(star: Star, level: i32) -> Meso {
    let level_factor: f32 = level.pow(3) as f32;
    let star_factor: f32 = ((star + 1) as f32).powf(2.7);
    let denom: f32 = match star {
        n if 10 <= n && n <= 14 => 400,
        n if 15 <= n && n <= 17 => 120,
        n if 18 <= n && n <= 19 => 110,
        n if 20 <= n && n <= 22 => 100,
        _ => panic!(),
    } as f32;
    let cost = 1000f32 + (level_factor as f32) * star_factor / denom;
    // let cost_approx = round(cost as i32, UNIT) / UNIT;
    println!(
        "at level {} it costs {} mesos at {} stars, star-dependent factor is {}",
        level,
        pp(cost as u64),
        star,
        star_factor / denom
    );
    (star_factor / denom / BASE_STAR_FACTOR).into()
}
