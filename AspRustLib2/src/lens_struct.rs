#[repr(C)]
pub struct Side
{
	pub r: f64,
	pub c: f64,
	pub k: f64,
	pub ad: f64,
	pub ae: f64,
	pub surf_type: i32 // 0 = plane, 1 = sphere with or wo conic, 2 = poly asphere
}

#[repr(C)]
pub struct Lens
{
	pub diameter: f64,
	//material: String,
	//lens_type: String,
	pub ap: f64,
	pub n: f64,
	pub wl: f64,
	pub ct: f64,
	pub side1: Side,
	pub side2: Side,
	pub efl: f64,
	pub bfl: f64
}