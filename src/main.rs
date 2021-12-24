#![feature(map_first_last)]

use std::time::SystemTime;

// use serde::{Deserialize, Serialize};

use rusqlite::NO_PARAMS;
use serde_rusqlite::*;

use maplestory_calculator::calc::*;
use maplestory_calculator::consts::*;

// TODO: also calculate # booms (independently)
// TODO: add level-dependent cost tables
// TODO: add options for star-catching, event discounts
// TODO: command-line options
// TODO: save outputs to a file
// TODO: plot outputs
// TODO: write README

fn main() {
    let t = SystemTime::now();
    let table = calculate3(160, false);
    println!("Finished in {} seconds", t.elapsed().unwrap().as_secs_f32());
    println!("Writing to file...");

    let table: Vec<(Star, Star, u64, f64)> =
        table.into_iter()
        .map(
            |((start, target), v)| {
                v.into_iter().map(|(c, p)| (start, target, c, p)).collect::<Vec<_>>()
            })
        .flatten()
        .collect();

    let connection = rusqlite::Connection::open_in_memory().unwrap();
    connection.execute("CREATE TABLE data (start INT, target INT, cost INT, prob REAL)", NO_PARAMS).unwrap();

    table.into_iter().for_each(
        |tuple| {
            connection.execute("INSERT INTO data (start, target, cost, prob) VALUES (?, ?, ?, ?)",
                               &to_params(&tuple).unwrap().to_slice()).unwrap();
        });

    connection.backup(rusqlite::DatabaseName::Main, "output.db", None).expect("failed db backup");
}
