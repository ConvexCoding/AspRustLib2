use super::ray_vector::{Ray, Vector3D};
use super::lens_struct::{Lens, Side};


pub fn trace_full_ray(ray: &Ray, lens: &Lens, refocus: f64) -> Ray
{
    // Trace ray from srf 0 to first lens surface. The axial distance here should be zero.
    let p2 = translate_to_surface(&ray.pvector, &ray.edir, &lens.side1, 0.0);
    let n2 = calc_slope(&p2, &lens.side1);
    let e2 = calc_dir_sines(&ray.edir, &n2, 1.0, lens.n); // after refraction

    // Trace to Surface 2 after refraction
    let p3 = translate_to_surface(&p2, &e2, &lens.side2, lens.ct);
    let n3 = calc_slope(
        &Vector3D {
            x: p3.x,
            y: p3.y,
            z: p3.z - lens.ct,
        },
        &lens.side2,
    ); // adjust z for CT of lens
    let e3 = calc_dir_sines(&e2, &n3, lens.n, 1.0);

    // transfer ray to image plane
    let p4 = translate_to_flat(&p3, &e3, lens.ct + lens.bfl + refocus);

    Ray {
        pvector: p4,
        edir: e3,
    }
}


pub fn translate_to_surface(p0: &Vector3D, e0: &Vector3D, side: &Side, plane: f64) -> Vector3D
{
    match side.surf_type
    {
        0 => translate_to_flat(p0, e0, plane),
        _ =>
        {
            let mut zest1 = calc_sag(p0.x, p0.y, &side, 0.001);
            let mut u = (zest1 - p0.z) / e0.z;
            let mut p1 = p0.clone();
            let mut p2 = p0 + e0 * u;

            for _ in 0..10
            {
                if (&p1 - &p2).length() > 1e-4
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

            p2
        }
    }
}

pub fn translate_to_flat(p: &Vector3D, e: &Vector3D, zplane: f64) -> Vector3D
{
    let u = (zplane - p.z) / e.z;
    p + &(e * u)
}

pub fn calc_slope(p: &Vector3D, s: &Side) -> Vector3D
{
    let r = p.x * p.x + p.y * p.y;
    let q0 = p.z - s.ad * r * r - s.ae * r * r * r;
    let q1 = -4.0 * s.ad * r - 6.0 * s.ae * r * r;

    let dx = p.x * (-s.c - s.c * (s.k + 1.0) * q1 * q0 + q1);
    let dy = p.y * (-s.c - s.c * (s.k + 1.0) * q1 * q0 + q1);
    let dz = 1.0 - s.c * (s.k + 1.0) * q0;

    let n = Vector3D {
        x: dx,
        y: dy,
        z: dz,
    };
    //let f = -(s.c / 2.0) * r - (s.c / 2.0) * (s.k + 1.0) * q0 * q0 + q0;
    &n / n.length()
}

pub fn calc_sag(x: f64, y: f64, side: &Side, rtolforzero: f64) -> f64
{
    let c = if side.r.abs() > rtolforzero
    {
        1.0 / side.r
    }
    else
    {
        0.0
    };

    let r2 = x * x + y * y;
    let sqrtvalue = 1.0 - (1.0 + side.k) * c * c * r2;

    if sqrtvalue < 0.0
    {
        0.0
    }
    else
    {
        c * r2 / (1.0 + sqrtvalue.sqrt()) + side.ad * r2 * r2 + side.ae * r2 * r2 * r2
    }
}

pub fn calc_dir_sines(ein: &Vector3D, ndir: &Vector3D, nin: f64, nout: f64) -> Vector3D
{
    let alpha = ein.dot_product(ndir);
    let a = 1.0;
    let b = 2.0 * alpha;
    let c = 1.0 - (nout * nout) / (nin * nin);
    let sol2 = (-b + (b * b - 4.0 * a * c).sqrt()) / (2.0 * a);
    let ep = ein + ndir * sol2;
    &ep / ep.length()
}
