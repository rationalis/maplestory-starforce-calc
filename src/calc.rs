use crate::consts::*;
use crate::distr::*;

use rustc_hash::FxHashMap;

pub fn calculate3(level: i32) {
    let cost = &(*COST); // just to trigger printlns
    let mut table: FxHashMap<(Star, Star), Distr> = FxHashMap::default();
    table.insert((12, 12), Distr::zero());
    let update = |table: &mut FxHashMap<_, _>, start: Star, target: Star, dist: Distr| {
        println!("{} -> {} dist size: {} expected cost: {}", start, target, dist.dist.len(), dist.expected_cost());
        // if start == 21 {
        //     dbg!(&dist.dist[0..100]);
        // }
        table.insert((start, target), dist);
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

            if down > 0. {
                let cost_below = COST[(start - 1) as usize];
                let [_, _, downdown, downboom] = PROBS_F64[(start - 11) as usize];
                // if we can go down then we need the number of failures
                if downdown == 0. {
                    let mut fails = base.clone();
                    fails.shift(-1);
                    let ds_cost = table.get(&(start-1, start)).unwrap();
                    dist = dist.add(&fails.product(&ds_cost));
                } else {
                    let cost_two_below = COST[(start - 2) as usize];
                    let [downups, downstays, downdowns, downbooms] = Distr::downs(start);
                    let downups = slice_to_distr(&downups);
                    let downstays = slice_to_distr(&downstays);
                    let downdowns = slice_to_distr(&downdowns);
                    let downbooms = slice_to_distr(&downbooms);
                    let mut du_cost = downups;
                    du_cost.scale(cost_below);
                    dist = dist.add(&du_cost);
                    let mut ds_cost = table.get(&(start - 1, start)).unwrap().clone();
                    ds_cost.shift(cost_below);
                    ds_cost = ds_cost.product(&downstays);
                    dist = dist.add(&ds_cost);
                    let mut dd_cost = table.get(&(start - 1, start)).unwrap().clone();
                    dd_cost.shift(cost_below).shift(cost_two_below);
                    dd_cost = dd_cost.product(&downdowns);
                    dist = dist.add(&dd_cost);
                    if downboom > 0. {
                        let mut db_cost = table.get(&(12, start)).unwrap().clone();
                        db_cost.shift(cost_below);
                        db_cost = db_cost.product(&downbooms);
                        dist = dist.add(&db_cost);
                    }
                }
            }
            if boom > 0. {
                let booms = Distr::booms(up, boom);
                let boom_cost = table.get(&(12, start)).unwrap();
                dist = dist.add(&booms.product(boom_cost));
            }
            update(&mut table, start, target, dist);
        }
    }
}
