use std::ops::{Add, Mul, Sub};

use noisy_float::prelude::*;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Meso {
    pub repr: R32
}

impl Add for Meso {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            repr: self.repr + other.repr
        }
    }
}

impl Mul for Meso {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            repr: self.repr * other.repr
        }
    }
}

impl Sub for Meso {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            repr: self.repr - other.repr
        }
    }
}

impl From<f32> for Meso {
    fn from(f: f32) -> Meso {
        Meso {
            repr: r32(f)
        }
    }
}

impl From<f64> for Meso {
    fn from(f: f64) -> Meso {
        Meso {
            repr: r32(f as f32)
        }
    }
}

impl From<Meso> for f32 {
    fn from(m: Meso) -> f32 {
        m.repr.raw()
    }
}

impl From<Meso> for f64 {
    fn from(m: Meso) -> f64 {
        m.repr.raw() as f64
    }
}

