use std::ops;

#[derive(Clone)]
#[repr(C)]
pub struct Ray
{
    pub pvector: Vector3D,
    pub edir: Vector3D,
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
#[derive(Debug, Clone)]
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
