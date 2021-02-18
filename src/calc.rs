use crate::consts::*;
use crate::distr::*;

use rustc_hash::FxHashMap;

pub fn pp(num: u64) -> String {
    use format_num::format_num;
    format_num!("0.3s", num as f64)
}

pub fn calculate3(level: i32, safeguard: bool) {
    let level = LEVELS.iter().position(|&e| e == level).unwrap();
    lazy_static::initialize(&BIN_SUMS);
    lazy_static::initialize(&COST);
    let mut table: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    let mut table_chance_time: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    table.insert((12, 12), Distr::zero());
    let update = |table: &mut FxHashMap<_, _>, start: Star, target: Star, dist: Distr| {
        let s = dist.stats();
        println!("{} -> {} dist size: {} mean: {} stddev: {}", start, target, dist.dist.len(), pp(s.0), pp(s.1));
        let q = dist.quartiles();
        println!("quartiles: {} {} {} {} {}", pp(q.0), pp(q.1), pp(q.2), pp(q.3), pp(q.4));
        table.insert((start, target), dist);
    };
    let dist_below = |table: &FxHashMap<_, Distr>, table_chance: &FxHashMap<_, _>, start: Star| {
        let key = &(start-1, start);
        table_chance.get(key)
            .or(table.get(key))
            .unwrap()
            .clone()
    };
    for target in 11..STAR_LIMIT+1 {
        for start in (10..target).rev() {
            if start != target - 1 {
                let dist1 = table.get(&(start, start+1)).unwrap();
                let dist2 = table.get(&(start+1, target)).unwrap();
                let dist = dist1.add(dist2);
                update(&mut table, start, target, dist);
                continue;
            }
            let cost = COST[level][start as usize];
            let [up, stay, mut down, mut boom] = PROBS_F64[(start - 10) as usize];
            if down == 0. && boom == 0. {
                let mut dist = Distr::geom(up);
                dist.scale(COST[level][start as usize]);
                update(&mut table, start, target, dist);
                continue;
            }

            if safeguard && start <= 16 {
                down += boom;
                boom = 0.0;
            }

            let joint = Distr::sim([up, stay, down, boom]);

            let cost_below = COST[level][(start - 1) as usize];
            let mut dist = PartialDistr::default();
            let mut dist_chance_time = PartialDistr::default();
            let mut keys: Vec<_> = joint.keys().collect();
            keys.sort_unstable();
            let keys = keys;
            let mut downs_dists = Vec::new();
            let mut booms_dists = Vec::new();
            downs_dists.push(Distr::zero());
            booms_dists.push(Distr::zero());
            let base_downs_dist = dist_below(&table, &table_chance_time, start);
            let base_booms_dist = if boom > 0.0 {
                table.get(&(12, start)).unwrap().clone()
            } else {
                Distr::zero()
            };
            downs_dists.push(base_downs_dist);
            if boom > 0.0 {
                booms_dists.push(base_booms_dist);
            }
            for key in keys {
                let &(downs, booms) = key;
                let (downs, booms) = (downs as usize, booms as usize);
                if downs == downs_dists.len() {
                    downs_dists.push(downs_dists.last().unwrap() + &downs_dists[1]);
                }
                if booms == booms_dists.len() {
                    booms_dists.push(booms_dists.last().unwrap() + &booms_dists[1]);
                }
                let booms_plus_downs = &downs_dists[downs] + &booms_dists[booms];
                let stays = joint.get(key).unwrap() * cost;
                if down > 0.0 && start < 21 {
                    let bpd_chance_time = if downs > 0 {
                        &(&downs_dists[downs - 1] + cost_below) + &booms_dists[booms]
                    } else {
                        booms_dists[booms].clone()
                    };
                    let total_cost_chance_time = bpd_chance_time + &stays;
                    dist_chance_time.mix(total_cost_chance_time);
                }
                let total_cost = booms_plus_downs + &stays;
                dist.mix(total_cost);
            }

            if start < 21 && down > 0.0 {
                update(&mut table_chance_time, start, target, dist_chance_time.into());
            }
            update(&mut table, start, target, dist.into());
        }
    }
}
