use crate::consts::*;
use crate::distr::*;

use rustc_hash::FxHashMap;

pub fn calculate3(level: i32) -> FxHashMap<(Star, Star), Vec<(u64, f64)>> {
    let cost = &(*COST); // just to trigger printlns
    let mut table: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    let mut table_chance_time: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    table.insert((12, 12), Distr::zero());
    let update = |table: &mut FxHashMap<_, _>, start: Star, target: Star, dist: Distr| {
        println!("{} -> {} dist size: {} expected cost: {}", start, target, dist.dist.len(), dist.expected_cost());
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
                let dist1 = &table[&(start, start+1)];
                let dist2 = &table[&(start+1, target)];
                let dist = dist1.add(dist2);
                update(&mut table, start, target, dist);
                continue;
            }
            let cost = COST[start as usize];
            let [up, _stay, down, boom] = PROBS_F64[(start - 10) as usize];
            if down == 0. && boom == 0. {
                let mut dist = Distr::geom(up);
                dist.scale(COST[start as usize]);
                update(&mut table, start, target, dist);
                continue;
            }

            let joint = Distr::sim(start);

            let cost_below = COST[(start - 1) as usize];
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
                table[&(12, start)].clone()
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
                print!("chance time: ");
                update(&mut table_chance_time, start, target, dist_chance_time.into());
            }
            update(&mut table, start, target, dist.into());
        }
    }
    table.into_iter()
        .map(|(k, v)|
             (k, Into::<Vec<(u64, f64)>>::into(v)))
        .collect()
}
