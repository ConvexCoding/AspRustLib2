use std::ops;

pub const CPROPV: Vector3D = Vector3D {
    x: 0.0,
    y: 0.0,
    z: 1.0,
};

#[derive(Clone)]
#[repr(C)]
pub struct Ray
{
    pub pvector: Vector3D,
    pub edir: Vector3D,
}

impl Ray
{
    pub fn calc_aoi_lsa(&self) -> (f64, f64)
    {
        let aoi = self.edir.dot_product(&CPROPV).acos();
        let lsa = -1.0 * (self.pvector.x.powi(2) + self.pvector.y.powi(2)).sqrt() / aoi.tan();
        (aoi, lsa)
    }
}

#[derive(Clone)]
#[repr(C)]
pub struct WFE_Ray
{
    pub rstart: Ray,
    pub rend: Ray,
    pub opd: f64,
    pub lsa: f64,
    pub ix: i32,
    pub iy: i32,
    pub isvalid: bool,
}

#[derive(Copy, Clone)]
#[repr(C)]
pub struct WFE_Stats
{
    pub minopd: f64,
    pub maxopd: f64,
    pub varirms: f64,
}

#[repr(C)]
#[derive(Debug, Clone, PartialEq)]
pub struct Vector3D
{
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Vector3D
{
    pub fn length(&self) -> f64
    {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
    }

    pub fn _lengthsquared(&self) -> f64
    {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    pub fn dot_product(&self, other: &Vector3D) -> f64
    {
        self.x * other.x + self.y * other.y + self.z * other.z
    }
}

// use this crates macro to more easily implement all the operator overloads.
// https://docs.rs/impl_ops/latest/impl_ops/
// we are always creating a new vector as result of the operation and need to handle the ref cases

impl_op!(+|a: Vector3D, b: Vector3D| -> Vector3D { &a + &b });
impl_op!(+|a: &Vector3D, b: Vector3D| -> Vector3D { a + &b });
impl_op!(+|a: Vector3D, b: &Vector3D| -> Vector3D { &a + b });
impl_op!(+|a: &Vector3D, b: &Vector3D| -> Vector3D {
    Vector3D {
        x: a.x + b.x,
        y: a.y + b.y,
        z: a.z + b.z,
    }
});

impl_op!(-|a: Vector3D, b: Vector3D| -> Vector3D { &a - &b });
impl_op!(-|a: &Vector3D, b: Vector3D| -> Vector3D { a - &b });
impl_op!(-|a: Vector3D, b: &Vector3D| -> Vector3D { &a - b });
impl_op!(-|a: &Vector3D, b: &Vector3D| -> Vector3D {
    Vector3D {
        x: a.x - b.x,
        y: a.y - b.y,
        z: a.z - b.z,
    }
});

// we dont actually use this, maybe * should be dot product?
// impl_op!(*|a: Vector3D, b: Vector3D| -> Vector3D { &a * &b });
// impl_op!(*|a: &Vector3D, b: Vector3D| -> Vector3D { a * &b });
// impl_op!(*|a: Vector3D, b: &Vector3D| -> Vector3D { &a * b });
// impl_op!(*|a: &Vector3D, b: &Vector3D| -> Vector3D {
//     Vector3D {
//         x: a.x * b.x,
//         y: a.y * b.y,
//         z: a.z * b.z,
//     }
// });

impl_op!(*|a: Vector3D, b: f64| -> Vector3D { &a * b });
impl_op!(*|a: &Vector3D, b: f64| -> Vector3D {
    Vector3D {
        x: a.x * b,
        y: a.y * b,
        z: a.z * b,
    }
});

impl_op!(/|a: Vector3D, b: f64| -> Vector3D { &a / b });
impl_op!(/|a: &Vector3D, b: f64| -> Vector3D {
    Vector3D {
        x: a.x / b,
        y: a.y / b,
        z: a.z / b,
    }
});

#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn add_vectors()
    {
        let a = Vector3D {
            x: 1.0,
            y: 1.0,
            z: 1.0,
        };
        let b = Vector3D {
            x: 2.0,
            y: 2.0,
            z: 2.0,
        };

        assert_eq!(
            a + b,
            Vector3D {
                x: 3.0,
                y: 3.0,
                z: 3.0,
            }
        )
    }

    #[test]
    fn sub_vectors()
    {
        let a = Vector3D {
            x: 3.0,
            y: 3.0,
            z: 3.0,
        };
        let b = Vector3D {
            x: 2.0,
            y: 2.0,
            z: 2.0,
        };

        assert_eq!(
            a - b,
            Vector3D {
                x: 1.0,
                y: 1.0,
                z: 1.0,
            }
        )
    }

    #[test]
    fn mul_vectors()
    {
        let a = Vector3D {
            x: 3.0,
            y: 3.0,
            z: 3.0,
        };

        assert_eq!(
            a * 3.0,
            Vector3D {
                x: 9.0,
                y: 9.0,
                z: 9.0,
            }
        )
    }

    #[test]
    fn div_vectors()
    {
        let a = Vector3D {
            x: 3.0,
            y: 3.0,
            z: 3.0,
        };

        assert_eq!(
            a / 2.0,
            Vector3D {
                x: 1.5,
                y: 1.5,
                z: 1.5,
            }
        )
    }

    #[test]
    fn dot_vectors()
    {
        let a = Vector3D {
            x: 1.0,
            y: 1.0,
            z: 1.0,
        };
        let b = Vector3D {
            x: 2.0,
            y: 2.0,
            z: 2.0,
        };

        assert_eq!(a.dot_product(&b), 6.0)
    }

    #[test]
    fn vector_length()
    {
        let a = Vector3D {
            x: 1.0,
            y: 1.0,
            z: 1.0,
        };
        assert_eq!(format!("{:.4}", a.length()), "1.7321")
    }
}
