#![feature(iter_map_while)]
#![feature(map_first_last)]
#![feature(option_result_unwrap_unchecked)]

use std::time::SystemTime;

use serde::{Deserialize, Serialize};

use maplestory_calculator::calc::*;

// TODO: also calculate # booms (independently)
// TODO: add level-dependent cost tables
// TODO: add options for star-catching, event discounts
// TODO: command-line options
// TODO: save outputs to a file
// TODO: plot outputs
// TODO: write README

fn main() {
    let t = SystemTime::now();
    calculate3(160);
    println!("Finished in {} seconds", t.elapsed().unwrap().as_secs_f32());
}
