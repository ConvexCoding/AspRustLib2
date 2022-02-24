use std::ops::{Add, Sub, Mul, Div};

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct Vector3D 
{
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[derive(Copy, Clone)]
#[repr(C)]
pub struct Ray
{
	pub pvector: Vector3D,
	pub edir:    Vector3D,
}

#[derive(Copy, Clone)]
#[repr(C)]
pub struct WFE_Ray
{

	pub rstart: Ray,
	pub rend: Ray,
	pub opd: f64,
	pub lsa: f64,
	pub ix: i32,
	pub iy: i32,
	pub isvalid: i32
}

#[derive(Copy, Clone)]
#[repr(C)]
pub struct WFE_Stats
{
	pub minopd: f64,
	pub maxopd: f64,
	pub varirms: f64
}

impl Vector3D
{
	pub fn length(self) -> f64
	{
		(self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
	}

	pub fn _lengthsquared(self) -> f64
	{
		self.x * self.x + self.y * self.y + self.z * self.z
	}
}

impl Add for Vector3D
{
	type Output = Self;
    fn add(self, other: Self) -> Self 
	{
        Self {x: self.x + other.x, y: self.y + other.y, z: self.z + other.z}
	}
}

impl Sub for Vector3D 
{
	type Output = Self;
    fn sub(self, other: Self) -> Self 
	{
        Self {x: self.x - other.x, y: self.y - other.y, z: self.z - other.z}
	}
}

impl Mul<f64> for Vector3D
{
	type Output = Self;
    fn mul(self, other: f64) -> Self 
	{
        Self {x: self.x * other, y: self.y * other, z: self.z * other}
	}
}

impl Mul<Vector3D> for Vector3D
{
	type Output = Self;
    fn mul(self, other: Self) -> Self 
	{
        Self {x: self.x * other.x, y: self.y * other.y, z: self.z * other.z}
	}
}

impl Div<f64> for Vector3D
{
	type Output = Self;
    fn div(self, other: f64) -> Self 
	{
        Self {x: self.x / other, y: self.y / other, z: self.z / other}
	}
}

pub fn dot_product(v1: Vector3D, v2: Vector3D) -> f64
{
	let x = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	return x;
}