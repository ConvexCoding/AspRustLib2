use crate::{Ray, Vector3D, dot_product, Side};

pub const  CPROPV: Vector3D = Vector3D{x: 0f64, y: 0f64, z: 1.0f64};

pub fn translate_to_surface(p0: Vector3D, e0: Vector3D, side: &Side, plane: f64) -> Vector3D
{
	if side.surf_type == 0
	{
		let pprime = translate_to_flat(p0, e0, plane);
		return pprime;
	}

	let mut zest1 = calc_sag(p0.x, p0.y, &side, 0.001);
	let mut u = (zest1 - p0.z) / e0.z;
	let mut p1 = p0.clone();
	let mut p2 = p0 + e0 * u;

	for _i in 0..10
	{
		if (p1 - p2).length() > 1e-4f64
		{
			p1 = p2;
			zest1 = calc_sag(p1.x, p1.y, &side, 0.001) + plane;
			u = (zest1 - p0.z) / e0.z;
			p2 = p0 + e0 * u;
		}
		else
		{
			break;
		}
	}

	return p2;
}

pub fn translate_to_flat(p: Vector3D, e: Vector3D, zplane: f64) -> Vector3D
{
	let u = (zplane - p.z) / e.z;
	let pprime = p + e * u;
	return pprime;
}


pub fn calc_aoi_lsa(ray: Ray) -> (f64, f64)
{
	//Vector3D{x: 0f64, y: 0f64, z: 1f64}
	let aoi = dot_product(ray.edir, CPROPV).acos();
	let lsa = -1.0f64 * (ray.pvector.x * ray.pvector.x + ray.pvector.y * ray.pvector.y).sqrt() / aoi.tan();
	return (aoi, lsa);
}

pub fn calc_slope(p: Vector3D, s: &Side) -> Vector3D 
{
	let r = p.x * p.x + p.y * p.y;
	let q0 = p.z - s.ad * r * r - s.ae * r * r * r;
	let q1 = -4.0 * s.ad * r - 6.0 * s.ae * r * r;

	let dx = p.x * (-s.c - s.c * (s.k + 1.0) * q1 * q0 + q1);
	let dy = p.y * (-s.c - s.c * (s.k + 1.0) * q1 * q0 + q1);
	let dz = 1.0 - s.c * (s.k + 1.0) * q0;

	let mut n = Vector3D{x: dx, y: dy, z: dz};
	n = n / n.length();
	//let f = -(s.c / 2.0) * r - (s.c / 2.0) * (s.k + 1.0) * q0 * q0 + q0;
	return n;
}

pub fn calc_sag(x: f64, y: f64, side: &Side, rtolforzero: f64) -> f64
{
	
	let mut c = 0.0;
	if side.r.abs() > rtolforzero
	{
		c = 1.0 / side.r;
	}

	let r2 = x * x + y * y;
	let sqrtvalue = 1.0 - (1.0 + side.k) * c * c * r2;

	if sqrtvalue < 0.0
	{
		return 0.0;
	}
	else
	{
		return c * r2 / (1.0 + sqrtvalue.sqrt()) + side.ad * r2 * r2 + side.ae * r2 * r2 * r2;
	}
}


pub fn calc_dir_sines(ein: Vector3D, ndir: Vector3D, nin: f64, nout: f64) -> Vector3D
{
	let alpha = dot_product(ein, ndir);
	let a = 1.0;
	let b = 2.0 * alpha;
	let c = 1.0 - (nout * nout) / (nin * nin);
	let sol2 = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
	let mut ep = ein + ndir * sol2;
	ep = ep / ep.length();
	return ep;
}