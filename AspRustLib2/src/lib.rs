use rand::Rng;

use crate::lens_struct::{Lens, Side};
mod lens_struct;

use crate::ray_vector::{Ray, Vector3D, WFE_Ray, WFE_Stats};
mod ray_vector;

use crate::rtm_subs::{
    calc_aoi_lsa, calc_dir_sines, calc_slope, translate_to_flat, translate_to_surface,
};
mod rtm_subs;

//extern crate regex;
//use regex::Regex;

#[repr(C)]
pub struct GenRays
{
    baserays: u32,
    num_angles: u32,
    half_ap: f64,
    half_ang: f64,
}

pub struct OpdResults
{
    pub rfinal: Ray,
    pub opd: f64,
    pub lsa: f64,
}

pub const CPROPV: Vector3D = Vector3D {
    x: 0f64,
    y: 0f64,
    z: 1.0f64,
};

//
// ===================================================================================
// External Functions visible over FFI
//

// +++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++

// method used in calculating opd
// single pvector and edir Vector3D's
// in the constrcution of WFE maps - AspGens
// calls ray trace and calc_aoi_lsa

#[no_mangle]
pub extern "C" fn rcalc_wfe(p0: Vector3D, e0: Vector3D, lens: &Lens, refocus: f64) -> f64
{
    let sqr2 = (2.0_f64).sqrt();
    let p1 = Vector3D {
        x: (p0.x / sqr2),
        y: (p0.y / sqr2),
        z: 1.0f64,
    };

    let rsq = p0.x * p0.x + p0.y * p0.y;
    let rsqsq = rsq * rsq;

    let rm = trace_ray(p0, e0, lens, 0.0f64);
    //let ym = rm.pvector;
    let (ymaoi, ymlsa) = calc_aoi_lsa(rm);

    let rz = trace_ray(p1, e0, lens, 0.0f64);
    //let yz = rz.pvector;
    let (_yzaoi, yzlsa) = calc_aoi_lsa(rz);

    //let rfinal = trace_ray(p0, e0, lens, refocus);
    //let yfinal = rfinal.pvector;
    //let (yfaoi, yflsa) = calc_aoi_lsa(rfinal);

    let a = (4.0f64 * yzlsa - ymlsa) / rsq;
    let b = (2.0f64 * ymlsa - 4.0f64 * yzlsa) / rsqsq;

    1000.0f64
        * (ymaoi.sin() * ymaoi.sin() / 2.0f64)
        * (refocus - a * rsq / 2.0f64 - b * rsqsq / 3.0f64)
        / lens.wl
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++

// method used to generate pupil array and calculate
// stats on traced rays
// this method accepts the WFE_Ray array pointers and
// modifies or populates the data the arrys with
// start ray, end ray, lsa, array position and validity
// method calls ray trace and calc_aoi_lsa methods

#[no_mangle]
extern "C" fn gen_and_trace_wfe_rays(
    ptr: *mut WFE_Ray,
    npts: usize,
    loopsize: usize,
    lens: Lens,
    refocus: f64,
) -> WFE_Stats
{
    let din: &mut [WFE_Ray] = unsafe { std::slice::from_raw_parts_mut(ptr, npts) };

    gen_wfe_rays(lens.ap, npts, din, loopsize);

    let mut wstats = WFE_Stats {
        minopd: 1e20f64,
        maxopd: -1e20f64,
        varirms: 0f64,
    };

    let mut xsum = 0.0f64;
    let mut xsumsq = 0.0f64;
    let mut cts: f64 = 1.0;

    for i in 0..npts
    {
        if din[i].isvalid == 1
        {
            calc_wfe_ray(&mut din[i], &lens, refocus);
            if !din[i].opd.is_nan()
            {
                if din[i].opd > wstats.maxopd
                {
                    wstats.maxopd = din[i].opd;
                }
                if din[i].opd < wstats.minopd
                {
                    wstats.minopd = din[i].opd;
                }

                xsum += din[i].opd;
                xsumsq += din[i].opd * din[i].opd;
                cts += 1.0f64;
            }
        }
    }
    wstats.varirms = ((xsumsq - xsum * xsum / cts) / (cts - 1f64)).sqrt();

    wstats
}

fn gen_wfe_rays(apert: f64, _size: usize, din: &mut [WFE_Ray], loopsize: usize)
{
    let mut x: f64;
    let mut y: f64;
    let diag = apert * apert;

    let step = 2f64 * apert / (loopsize - 1) as f64;
    let mut count = 0usize;

    //let dirvector = Vector3D{x: 0.0f64, y: 0.0f64, z: 1.0f64};

    for row in 0..loopsize
    {
        y = apert - row as f64 * step;
        for col in 0..loopsize
        {
            x = -apert + col as f64 * step;

            din[count].rstart.pvector = Vector3D { x: x, y: y, z: 0.0 };
            din[count].rstart.edir = CPROPV;
            din[count].ix = col as i32;
            din[count].iy = row as i32;

            if diag > x * x + y * y
            {
                din[count].isvalid = 1i32;
            }
            else
            {
                din[count].isvalid = 0i32;
            }
            count = count + 1
        }
    }
}

fn calc_wfe_ray(wferay: &mut WFE_Ray, lens: &Lens, refocus: f64)
{
    let sqr2 = (2.0_f64).sqrt();

    let p0 = wferay.rstart.pvector;
    let e0 = wferay.rstart.edir;

    let p1 = Vector3D {
        x: (p0.x / sqr2),
        y: (p0.y / sqr2),
        z: 1.0f64,
    };

    let rsq = p0.x * p0.x + p0.y * p0.y;

    if rsq < 0.0000001f64
    {
        //wferay.rend = Ray{pvector: Vector3D{x: 0f64, y: 0f64, z: lens.ct + lens.bfl + refocus}, edir: Vector3D{x: 0f64, y: 0f64, z: 1f64}};
        wferay.rend = Ray {
            pvector: Vector3D {
                x: 0f64,
                y: 0f64,
                z: lens.ct + lens.bfl + refocus,
            },
            edir: CPROPV,
        };
        wferay.lsa = 0f64;
        wferay.opd = 0f64;
    }

    let rsqsq = rsq * rsq;

    let rm = trace_ray(p0, e0, lens, 0.0f64);
    //let ym = rm.pvector;
    let (ymaoi, ymlsa) = calc_aoi_lsa(rm);

    let rz = trace_ray(p1, e0, lens, 0.0f64);
    //let yz = rz.pvector;
    let (_yzaoi, yzlsa) = calc_aoi_lsa(rz);

    let rfinal = trace_ray(p0, e0, lens, refocus);
    //let yfinal = rfinal.pvector;
    let (_yfaoi, yflsa) = calc_aoi_lsa(rfinal);
    wferay.rend = rfinal;
    wferay.lsa = yflsa;

    let a = (4.0f64 * yzlsa - ymlsa) / rsq;
    let b = (2.0f64 * ymlsa - 4.0f64 * yzlsa) / rsqsq;

    wferay.opd = 1000.0f64
        * (ymaoi.sin() * ymaoi.sin() / 2.0f64)
        * (refocus - a * rsq / 2.0f64 - b * rsqsq / 3.0f64)
        / lens.wl;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++

#[no_mangle]
pub extern "C" fn tracerays(
    ptrin: *mut Ray,
    prtout: *mut Ray,
    npts: usize,
    lens: Lens,
    refocus: f64,
    aposi: usize,
) -> f64
{
    let din: &mut [Ray] = unsafe { std::slice::from_raw_parts_mut(ptrin, npts) };
    let dout: &mut [Ray] = unsafe { std::slice::from_raw_parts_mut(prtout, npts) };

    for i in 0..npts
    {
        let tr = trace_ray(din[i].pvector, din[i].edir, &lens, refocus);
        dout[i].pvector = tr.pvector;
        dout[i].edir = tr.edir;
    }

    dout[aposi].pvector.x
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++

#[no_mangle]
pub extern "C" fn gen_trace_rays(
    gr: GenRays,
    ptrin: *mut Ray,
    prtout: *mut Ray,
    npts: usize,
    lens: Lens,
    refocus: f64,
) -> bool
{
    let din: &mut [Ray] = unsafe { std::slice::from_raw_parts_mut(ptrin, npts) };
    let dout: &mut [Ray] = unsafe { std::slice::from_raw_parts_mut(prtout, npts) };

    gen_random_rays(gr, din);

    for i in 0..npts
    {
        let tr = trace_ray(din[i].pvector, din[i].edir, &lens, refocus);
        dout[i].pvector = tr.pvector;
        dout[i].edir = tr.edir;
    }

    true
}

#[no_mangle]
pub fn gen_random_rays(gr: GenRays, din: &mut [Ray])
{
    let mut x: f64;
    let mut y: f64;
    let mut xdir: f64;
    let mut ydir: f64;

    let diag = gr.half_ap * gr.half_ap;
    let anglediag = gr.half_ang * gr.half_ang;

    let mut count: usize = 0;

    let mut rng = rand::thread_rng();

    for _ in 0..gr.baserays
    {
        x = rng.gen_range(-gr.half_ap, gr.half_ap);
        y = rng.gen_range(-gr.half_ap, gr.half_ap);
        while (x * x + y * y) > diag
        {
            x = rng.gen_range(-gr.half_ap, gr.half_ap);
            y = rng.gen_range(-gr.half_ap, gr.half_ap);
        }
        let pvbase = Vector3D { x: x, y: y, z: 0.0 };

        for _j in 0..gr.num_angles
        {
            xdir = rng.gen_range(-gr.half_ang, gr.half_ang);
            ydir = rng.gen_range(-gr.half_ang, gr.half_ang);
            while (xdir * xdir + ydir * ydir) > anglediag
            {
                xdir = rng.gen_range(-gr.half_ang, gr.half_ang);
                ydir = rng.gen_range(-gr.half_ang, gr.half_ang);
            }
            let edir = Vector3D {
                x: xdir,
                y: ydir,
                z: (1.0 - xdir * xdir - ydir * ydir).sqrt(),
            };
            din[count].pvector = pvbase;
            din[count].edir = edir;
            count += 1;
        }
    }
}

#[no_mangle]
pub extern "C" fn trace_ray(p0: Vector3D, e0: Vector3D, lens: &Lens, refocus: f64) -> Ray
{
    // Trace ray from srf 0 to first lens surface. The axial distance here should be zero.
    let p2 = translate_to_surface(p0, e0, &lens.side1, 0.0);
    let n2 = calc_slope(p2, &lens.side1);
    let e2 = calc_dir_sines(e0, n2, 1.0, lens.n); // after refraction

    // Trace to Surface 2 after refraction
    let p3 = translate_to_surface(p2, e2, &lens.side2, lens.ct);
    let n3 = calc_slope(
        Vector3D {
            x: p3.x,
            y: p3.y,
            z: p3.z - lens.ct,
        },
        &lens.side2,
    ); // adjust z for CT of lens
    let e3 = calc_dir_sines(e2, n3, lens.n, 1.0);

    // transfer ray to image plane
    let p4 = translate_to_flat(p3, e3, lens.ct + lens.bfl + refocus);

    Ray {
        pvector: p4,
        edir: e3,
    }
}

#[no_mangle]
pub extern "C" fn trace_full_ray(ray: Ray, lens: &Lens, refocus: f64) -> Ray
{
    // Trace ray from srf 0 to first lens surface. The axial distance here should be zero.
    let p2 = translate_to_surface(ray.pvector, ray.edir, &lens.side1, 0.0);
    let n2 = calc_slope(p2, &lens.side1);
    let e2 = calc_dir_sines(ray.edir, n2, 1.0, lens.n); // after refraction

    // Trace to Surface 2 after refraction
    let p3 = translate_to_surface(p2, e2, &lens.side2, lens.ct);
    let n3 = calc_slope(
        Vector3D {
            x: p3.x,
            y: p3.y,
            z: p3.z - lens.ct,
        },
        &lens.side2,
    ); // adjust z for CT of lens
    let e3 = calc_dir_sines(e2, n3, lens.n, 1.0);

    // transfer ray to image plane
    let p4 = translate_to_flat(p3, e3, lens.ct + lens.bfl + refocus);

    Ray {
        pvector: p4,
        edir: e3,
    }
}
