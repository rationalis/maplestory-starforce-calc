use crate::consts::*;
use crate::distr;
use crate::distr::*;

use std::borrow::Cow;
use std::thread;
use std::sync::mpsc;

use rustc_hash::FxHashMap;

pub fn calculate3(level: i32, safeguard: bool) -> Vec<((Star, Star), Vec<(u64, f64)>)> {
    // TODO: add back reporting
    // let (to_reporter, reporter_receiver) = mpsc::channel::<(Star, Star, Distr)>();
    // let _reporter = thread::spawn(move || {
    //     while let Ok((start, target, dist)) = reporter_receiver.recv() {
    //         let s = dist.stats();
    //         println!(
    //             "{} -> {} dist size: {} mean: {} stddev: {}",
    //             start,
    //             target,
    //             dist.dist.len(),
    //             pp(s.0),
    //             pp(s.1)
    //         );
    //         let q = dist.quartiles();
    //         println!(
    //             "quartiles: {} {} {} {} {}",
    //             pp(q.0),
    //             pp(q.1),
    //             pp(q.2),
    //             pp(q.3),
    //             pp(q.4)
    //         );
    //     }
    // });

    // let level = LEVELS.iter().position(|&e| e == level).unwrap();
    lazy_static::initialize(&BIN_SUMS);
    lazy_static::initialize(&COST);
    // TODO: add separate tables for implicits -- this will save on FFTs
    let mut table: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    let mut table_chance_time: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    // for booming
    table.insert((12, 12), Distr::zero());
    let update = |table: &mut FxHashMap<_, _>, start: Star, target: Star, dist: Distr| {
        // to_reporter.send((start, target, dist.clone())).unwrap();
        table.insert((start, target), dist);
    };
    let dist_below = |table: &FxHashMap<_, Distr>, table_chance: &FxHashMap<_, _>, start: Star| {
        let key = &(start - 1, start);
        table_chance.get(key).or(table.get(key)).unwrap().clone()
    };
    for target in 11..STAR_LIMIT + 1 {
        for start in (10..target).rev() {
            if start != target - 1 {
                let dist1 = &table[&(start, start + 1)];
                let dist2 = &table[&(start + 1, target)];
                // TODO: with a separate implicits table, we can replace this with a normal conv
                let dist = (&dist1.characteristic() + &dist2.characteristic()).explicit();
                update(&mut table, start, target, dist);
                continue;
            }
            let cost: Meso = COST[start as usize];

            let [up, stay, mut down, mut boom] = PROBS_F64[(start - 10) as usize];
            if down == 0. && (boom == 0. || safeguard) {
                let mut dist = distr::geom(up);
                dist *= f(COST[start as usize]);
                update(&mut table, start, target, dist);
                continue;
            }

            if safeguard && start <= 16 {
                down += boom;
                boom = 0.0;
            }

            // un-mut
            let [up, stay, down, boom] = [up, stay, down, boom];

            let joint = distr::sim([up, stay, down, boom]);

            // let cost_below = COST[level][(start - 1) as usize];
            let cost_below = COST[(start - 1) as usize];
            let mut dist: ImplicitDistr = distr::all_empty();
            let mut dist_chance_time: ImplicitDistr = distr::all_empty();
            let mut keys: Vec<_> = joint.keys().collect();
            keys.sort_unstable();
            let keys = keys;
            let mut downs_dists: Vec<CharacteristicDistr> = Vec::new();
            let mut booms_dists: Vec<CharacteristicDistr> = Vec::new();
            // zeros (should NEVER actually be fft'd/convolved)
            downs_dists.push(Default::default());
            booms_dists.push(Default::default());
            // distributions of cost for single down/boom
            downs_dists.push(dist_below(&table, &table_chance_time, start).characteristic());
            if boom > 0.0 {
                booms_dists.push(table[&(12, start)].clone().characteristic());
            }
            // TODO: use max_downs, max_booms and separately precompute instead (parallel possible)
            // TODO: investigate -- computing directly from (downs, booms) pair sums to halve convolutions
            for key in keys {
                let &(downs, booms) = key;
                let (downs, booms) = (downs as usize, booms as usize);
                if downs == downs_dists.len() {
                    downs_dists.push(downs_dists.last().unwrap() + &downs_dists[1]);
                }
                if booms == booms_dists.len() {
                    booms_dists.push(booms_dists.last().unwrap() + &booms_dists[1]);
                }
                // avoid convolution on zero distr
                let booms_plus_downs  = match (downs, booms) {
                    (0, _) => Cow::Borrowed(&booms_dists[booms]),
                    (_, 0) => Cow::Borrowed(&downs_dists[downs]),
                    _ => Cow::Owned(&downs_dists[downs] + &booms_dists[booms])
                };

                let mut stays = joint.get(key).unwrap().clone();
                stays *= f(cost);
                let stays = &stays.characteristic();

                if down > 0.0 && start < 21 {
                    let bpd_chance_time = if downs > 0 {
                        Cow::Owned(&downs_dists[downs - 1] + &booms_dists[booms])
                    } else {
                        Cow::Borrowed(&booms_dists[booms])
                    };
                    let total_cost_chance_time = &*bpd_chance_time + stays;
                    let total_cost_chance_time = total_cost_chance_time.ifft().shift(cost_below);
                    dist_chance_time += &total_cost_chance_time;
                }
                let total_cost = &*booms_plus_downs + stays;
                dist += &total_cost;
            }

            if start < 21 && down > 0.0 {
                update(
                    &mut table_chance_time,
                    start,
                    target,
                    dist_chance_time.into(),
                );
            }
            update(&mut table, start, target, dist.into());
        }
    }
    table.into_iter()
        .map(|(k, v)| (k, v.star_factor_to_mesos(level)))
        .collect()
}
