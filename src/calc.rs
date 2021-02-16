use crate::consts::*;
use crate::distr::*;

use rustc_hash::FxHashMap;

pub fn calculate3(level: i32) {
    let cost = &(*COST); // just to trigger printlns
    let mut table: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    let mut table_chance_time: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    table.insert((12, 12), Distr::zero());
    let update = |table: &mut FxHashMap<_, _>, start: Star, target: Star, dist: Distr| {
        println!("{} -> {} dist size: {} expected cost: {}", start, target, dist.dist.len(), dist.expected_cost());
        // if start == 21 {
        //     dbg!(&dist.dist[0..100]);
        // }
        table.insert((start, target), dist);
    };
    let dist_below = |table: &FxHashMap<_, Distr>, table_chance: &FxHashMap<_, _>, start: Star| {
        let mut ds_cost = table_chance.get(&(start-1, start));
        if ds_cost.is_some() {
            println!("{}->{}->{} found chance time", start-1, start, start+1);
        } else {
            ds_cost = table.get(&(start-1, start));
        }
        let ds_cost = ds_cost.unwrap();
        ds_cost.clone()
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
            let [up, _stay, down, boom] = PROBS_F64[(start - 10) as usize];
            if down == 0. && boom == 0. {
                let mut dist = Distr::geom(up);
                dist.scale(COST[start as usize]);
                update(&mut table, start, target, dist);
                continue;
            }
            let mut dist = Distr::zero();
            let base = Distr::geom(up);
            // cost of attempts at start->target
            let mut base2 = base.clone();
            base2.scale(COST[start as usize]);
            dist = dist.add(&base2);

            let mut dist_chance_time = None;

            if down > 0. {
                let cost_below = COST[(start - 1) as usize];
                let [_, _, downdown, _] = PROBS_F64[(start - 11) as usize];
                // if we can go down then we need the number of failures
                // if downdown == 0. {
                if start < 21 {
                    let mut fails = base.clone();
                    fails.shift(-1);
                    let normal_cost = dist_below(&table, &table_chance_time, start);
                    let chance_time_cost = fails.chance_time(&normal_cost, cost_below);
                    dist_chance_time = Some(dist.add(&chance_time_cost));
                }
                let mut fails = base.clone();
                fails.shift(-1);
                let ds_cost = dist_below(&table, &table_chance_time, start);
                dbg!(start, fails.dist[0], fails.dist.len(), ds_cost.dist.len());
                let ds_cost = fails.product(&ds_cost);
                dist = dist.add(&ds_cost);
                dbg!(dist.dist.len());


                // } else {
                // if downdown > 0. {
                //     let cost_two_below = COST[(start - 2) as usize];
                //     let dd = Distr::downdowns(start);
                //     // subtract the start-1 -> start cost in downdowns (chance time)
                //     let mut dd_normal_cost = table.get(&(start - 2, start - 1)).unwrap().clone();
                //     dd_normal_cost = dd_normal_cost.product(&dd);
                //     dd_normal_cost.scale(-1);
                //     dist = dist.add(&dd_normal_cost);
                //     let mut chance_time_cost = Distr::constant(cost_two_below);
                //     chance_time_cost = chance_time_cost.product(&dd);
                //     dist = dist.add(&chance_time_cost);
                //     // dbg!(&dd_normal_cost.dist[0..10], &chance_time_cost.dist[0..10]);
                //     dbg!(dd_normal_cost.dist.len(), chance_time_cost.dist.len(), dist.dist.len());
                // }
            }
            if boom > 0. {
                let booms = Distr::booms(up, boom);
                let boom_cost = table.get(&(12, start)).unwrap();
                dist = dist.add(&booms.product(boom_cost));
            }
            if let Some(dist) = dist_chance_time {
                update(&mut table_chance_time, start, target, dist);
            }
            update(&mut table, start, target, dist);
        }
    }
}
