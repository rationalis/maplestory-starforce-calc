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
                if start < 21 {
                    // if we are under 21* and the down chance is >0, then chance
                    // time may trigger when we fall to this star and fall again
                    let mut fails = base.clone();
                    fails.shift(-1);
                    // the normal cost when falling to the star below this one -
                    // note that this is also affected by chance time so we make
                    // a call to dist_below
                    let normal_cost = dist_below(&table, &table_chance_time, start);
                    let chance_time_cost = fails.chance_time(&normal_cost, cost_below);
                    dbg!(&chance_time_cost.dist[0..100]);
                    dist_chance_time = Some(dist.add(&chance_time_cost));
                }
                let mut fails = base.clone();
                fails.shift(-1);
                // dbg!(start, &fails.dist);
                let ds_cost = dist_below(&table, &table_chance_time, start);
                let ds_cost = fails.product(&ds_cost);
                dbg!(&ds_cost.dist[0..100]);
                dist = dist.add(&ds_cost);
            }
            if boom > 0. {
                let booms = Distr::booms(up, boom);
                let boom_cost = table.get(&(12, start)).unwrap();
                let boom_cost = &booms.product(boom_cost);
                dist = dist.add(boom_cost);
                dist_chance_time = dist_chance_time.map(|d| d.add(boom_cost));
            }
            if let Some(dist) = dist_chance_time {
                print!("chance time: ");
                update(&mut table_chance_time, start, target, dist);
            }
            update(&mut table, start, target, dist);
        }
    }
}
