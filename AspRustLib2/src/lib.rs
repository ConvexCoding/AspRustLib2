#[macro_use]
extern crate impl_ops;

use rand::Rng;

mod tools;
use tools::rtm_subs::{trace_full_ray, calc_dir_sines, calc_slope, translate_to_flat, translate_to_surface};
use tools::ffttools::{gen_zero_2d, get_complex_vec, slicecore, find_max};
use tools::fft2d::rustfft;
use tools::ray_vector::{Ray, Vector3D, WFE_Ray, WFE_Stats, CPROPV};
use tools::lens_struct::{Lens, OptVariables};
use tools::optimize::opti_lens;

use std::f64::consts::SQRT_2 as sqr2;
use std::f64::consts::PI;


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

#[repr(C)]
pub struct ErrorStat
{
    peak: f64,
    valley: f64,
    pv: f64,
    average: f64,
    rms: f64,
}

//
// ===================================================================================
// External Functions visible over FFI
//

// +++++++++++   GEN CalcErrorFunc  ++++++++++++++++++

#[no_mangle]
extern "C" fn calc_err_func(
    ptr1: *mut f64,
    ptr2: *mut Ray,
    rayct: usize,
    lens: Lens,
    refocus: f64,
) -> f64
{
    let ray_ys: &mut [f64] = unsafe { std::slice::from_raw_parts_mut(ptr1, rayct) };
    let rayout: &mut [Ray] = unsafe { std::slice::from_raw_parts_mut(ptr2, rayct) };

    let mut sumsum: f64 = 0.0;

    let mut _mlens: Lens = lens;

    for i in 0..rayct{
        let ray: Ray = Ray{pvector: Vector3D{x: 0.0, y: ray_ys[i] * lens.ap, z: 0.0}, edir: CPROPV};
        rayout[i] = trace_full_ray(&ray, &lens, refocus);
        sumsum += rayout[i].pvector.y * rayout[i].pvector.y;
    }
    //mlens.side2.r = 2000.0;
    return (sumsum / rayct as f64).sqrt();
}

#[no_mangle]
extern "C" fn optimize_lens(
    ptr1: *mut f64,
    rayct: usize,
    lens: Lens,
    optv: OptVariables
) -> Lens
{
    let ray_ys: &mut [f64] = unsafe { std::slice::from_raw_parts_mut(ptr1, rayct) };

    /*
    let mut tracker: i32 = 0;
    if optv.iscc1on == 1 { tracker = -1;}
    if optv.isad1on == 1 { tracker = -2;}
    if optv.isae1on == 1 { tracker = -3;}
    if optv.iscc2on == 1 { tracker = -4;}
    if optv.isad2on == 1 { tracker = -5;}
    if optv.isae2on == 1 { tracker = -6;}
    */

    let mlens = opti_lens(&ray_ys, lens, &optv);
    return mlens;
}

#[no_mangle]
extern "C" fn calc_meritfunction_stats(
    ptr1: *mut f64,
    rayct: usize,
    lens: Lens,
    refocus: f64,
) -> ErrorStat
{
    let ray_ys: &mut [f64] = unsafe { std::slice::from_raw_parts_mut(ptr1, rayct) };

    let mut sum: f64 = 0.0;
    let mut sumsum: f64 = 0.0;
    let mut peak: f64 = -1.0e20;
    let mut valley: f64 = 1.0e20;

    let mut _mlens: Lens = lens;

    for i in 0..rayct{
        let ray: Ray = Ray{pvector: Vector3D{x: 0.0, y: ray_ys[i] * lens.ap, z: 0.0}, edir: CPROPV};
        let rimage = trace_full_ray(&ray, &lens, refocus);
        if rimage.pvector.y > peak {
            peak = rimage.pvector.y
        }
        if rimage.pvector.y < valley {
            valley = rimage.pvector.y
        }
        sum += rimage.pvector.y;
        sumsum += rimage.pvector.y * rimage.pvector.y;
    }
    let estat = ErrorStat{peak, valley, pv: (peak - valley), average: sum / rayct as f64, rms: (sumsum / rayct as f64).sqrt()};
    return estat;
}

#[no_mangle]
extern "C" fn calc_wfe_stats(
    ptr1: *mut f64,
    rayct: usize,
    lens: Lens,
    refocus: f64,
) -> ErrorStat
{
    let ray_ys: &mut [f64] = unsafe { std::slice::from_raw_parts_mut(ptr1, rayct) };

    let mut sum: f64 = 0.0;
    let mut sumsum: f64 = 0.0;
    let mut peak: f64 = -1.0e20;
    let mut valley: f64 = 1.0e20;

    let mut _mlens: Lens = lens;

    for i in 0..rayct{
        let wfe = calc_opd_slim(Vector3D{x: 0.0, y: ray_ys[i] * lens.ap, z: 0.0}, CPROPV, &lens, refocus);
        if wfe > peak {
            peak = wfe
        }
        if wfe < valley {
            valley = wfe
        }
        sum += wfe;
        sumsum += wfe * wfe;
    }
    let estat = ErrorStat{peak, valley, pv: (peak - valley), average: sum / rayct as f64, rms: (sumsum / rayct as f64).sqrt()};
    return estat;
}
// +++++++++++++++  GEN CalcErrorFunc END ++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++++++

//  private unsafe static extern double[] gen_psf(UIntPtr loopsize, UIntPtr paddedsize, UIntPtr psfsize, LensWOStrings lens, double refocus);
//
//  external fn callable from c# that will generate a psf given lens, ap, and various grid sizes
//  most of the work is down in folder tools, files fft2d.rs and ffttools.rs
//  one exception is that the wavefront is first constructed using the slightly optimized calc_opd fn added to this file

#[no_mangle]
extern "C" fn gen_psf(
    ptr: *mut f64,
    loopsize: usize,
    totalsize: usize,
    psfgridsize: usize,
    lens: Lens,
    refocus: f64,
) -> f64
{
    let psfdata: &mut [f64] = unsafe { std::slice::from_raw_parts_mut(ptr, psfgridsize * psfgridsize) };

    //gen_wfe_rays(lens.ap, totalsize, din, loopsize);

    let mut amp = gen_zero_2d(loopsize);
    let mut mask = gen_zero_2d(loopsize);

    let mut x: f64;
    let mut y: f64;
    let diag = lens.ap * lens.ap;

    let step = 2.0 * lens.ap / (loopsize - 1) as f64;
    
    for row in 0..psfgridsize
    {
        y = lens.ap - row as f64 * step;
        for col in 0..psfgridsize
        {
            x = -lens.ap + col as f64 * step;
            if (x * x + y * y) < diag {
                let p0 = Vector3D { x, y, z: 0.0 };
                amp[row][col] =   2.0 * PI * calc_opd_slim(p0, CPROPV, &lens, refocus);  // multiply by 2pi to scale for fft
                mask[row][col] = 1.0;
            }
            else { 
                amp[row][col] = 0.0;
                mask[row][col] = 0.0;
            }
        }
    }

    let mut data = get_complex_vec(&amp, &mask, totalsize);
    let mut datadl = get_complex_vec(&mask, &mask, totalsize);
    let datafull = rustfft(&mut data, &mut datadl);
    let _datatoc = slicecore(datafull, psfgridsize);

    for row in 0..psfgridsize
    {
        for col in 0..psfgridsize
        {
            psfdata[row * psfgridsize + col] = _datatoc[row][col];
        }
    }


    return find_max(_datatoc);
}


// this opd calc is used specifically for FFT2D calculations and is optimized for only opd
fn calc_opd_slim(p0: Vector3D, e0: Vector3D, lens: &Lens, refocus: f64) -> f64
{
    //let sqr2 = 2_f64.sqrt();

    let p1 = Vector3D {
        x: (p0.x / sqr2),
        y: (p0.y / sqr2),
        z: 1.0,
    };

    let rsq = p0.x * p0.x + p0.y * p0.y;

    if rsq < 1.0e-10 {
        return 0.0_f64;
    }

    let rsqsq = rsq * rsq;

    let rm = trace_ray(&p0, &e0, lens, 0.0);
    //let ym = rm.pvector;
    let (ymaoi, ymlsa) = rm.calc_aoi_lsa();

    let rz = trace_ray(&p1, &e0, lens, 0.0);
    //let yz = rz.pvector;
    let (_yzaoi, yzlsa) = rz.calc_aoi_lsa();

    let rfinal = trace_ray(&p0, &e0, lens, refocus);
    //let yfinal = rfinal.pvector;
    let (_yfaoi, _yflsa) = rfinal.calc_aoi_lsa();

    let a = (4.0 * yzlsa - ymlsa) / rsq;
    let b = (2.0 * ymlsa - 4.0 * yzlsa) / rsqsq;

    1000.0 * (ymaoi.sin() * ymaoi.sin() / 2.0) * (refocus - a * rsq / 2.0 - b * rsqsq / 3.0) / lens.wl

}


// +++++++++++++++   GEN PSF  End  +++++++++++++++++++

// +++++++++++++++++++++++++++++++++++++++++++++++++++
// +++++++++++++++++++++++++++++++++++++++++++++++++++

// method used in calculating opd
// single pvector and edir Vector3D's
// in the constrcution of WFE maps - AspGens
// calls ray trace and calc_aoi_lsa

#[no_mangle]
pub extern "C" fn rcalc_wfe(p0: Vector3D, e0: Vector3D, lens: &Lens, refocus: f64) -> f64
{
    let p1 = Vector3D {
        x: (p0.x / sqr2),
        y: (p0.y / sqr2),
        z: 1.0,
    };

    let rsq = p0.x * p0.x + p0.y * p0.y;
    let rsqsq = rsq * rsq;

    let rm = trace_ray(&p0, &e0, lens, 0.0);
    //let ym = rm.pvector;
    let (ymaoi, ymlsa) = rm.calc_aoi_lsa();

    let rz = trace_ray(&p1, &e0, lens, 0.0);
    //let yz = rz.pvector;
    let (_yzaoi, yzlsa) = rz.calc_aoi_lsa();

    //let rfinal = trace_ray(p0, e0, lens, refocus);
    //let yfinal = rfinal.pvector;
    //let (yfaoi, yflsa) = calc_aoi_lsa(rfinal);

    let a = (4.0 * yzlsa - ymlsa) / rsq;
    let b = (2.0 * ymlsa - 4.0 * yzlsa) / rsqsq;

    1000.0 * (ymaoi.sin() * ymaoi.sin() / 2.0) * (refocus - a * rsq / 2.0 - b * rsqsq / 3.0)
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
        minopd: 1e20,
        maxopd: -1e20,
        varirms: 0.0,
    };

    let mut xsum = 0.0;
    let mut xsumsq = 0.0;
    let mut cts = 1.0;

    for i in 0..npts
    {
        if din[i].isvalid
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
                cts += 1.0;
            }
        }
    }
    wstats.varirms = ((xsumsq - xsum * xsum / cts) / (cts - 1.0)).sqrt();

    wstats
}

fn gen_wfe_rays(apert: f64, _size: usize, din: &mut [WFE_Ray], loopsize: usize)
{
    let mut x: f64;
    let mut y: f64;
    let diag = apert * apert;

    let step = 2.0 * apert / (loopsize - 1) as f64;
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
            din[count].isvalid = diag > x * x + y * y;

            count += 1
        }
    }
}

fn calc_wfe_ray(wferay: &mut WFE_Ray, lens: &Lens, refocus: f64)
{
    let p0 = &wferay.rstart.pvector;
    let e0 = &wferay.rstart.edir;

    let p1 = Vector3D {
        x: (p0.x / sqr2),
        y: (p0.y / sqr2),
        z: 1.0,
    };

    let rsq = p0.x * p0.x + p0.y * p0.y;

    if rsq < 0.0000001
    {
        //wferay.rend = Ray{pvector: Vector3D{x: 0f64, y: 0f64, z: lens.ct + lens.bfl + refocus}, edir: Vector3D{x: 0f64, y: 0f64, z: 1f64}};
        wferay.rend = Ray {
            pvector: Vector3D {
                x: 0.0,
                y: 0.0,
                z: lens.ct + lens.bfl + refocus,
            },
            edir: CPROPV,
        };
        wferay.lsa = 0.0;
        wferay.opd = 0.0;
    }

    let rsqsq = rsq * rsq;

    let rm = trace_ray(&p0, &e0, lens, 0.0);
    //let ym = rm.pvector;
    let (ymaoi, ymlsa) = rm.calc_aoi_lsa();

    let rz = trace_ray(&p1, &e0, lens, 0.0);
    //let yz = rz.pvector;
    let (_yzaoi, yzlsa) = rz.calc_aoi_lsa();

    let rfinal = trace_ray(&p0, &e0, lens, refocus);
    //let yfinal = rfinal.pvector;
    let (_yfaoi, yflsa) = rfinal.calc_aoi_lsa();
    wferay.rend = rfinal;
    wferay.lsa = yflsa;

    let a = (4.0 * yzlsa - ymlsa) / rsq;
    let b = (2.0 * ymlsa - 4.0 * yzlsa) / rsqsq;

    wferay.opd =
        1000.0 * (ymaoi.sin() * ymaoi.sin() / 2.0) * (refocus - a * rsq / 2.0 - b * rsqsq / 3.0)
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
        //let tr = trace_ray(&din[i].pvector, &din[i].edir, &lens, refocus);
        let tr = trace_full_ray(&din[i], &lens, refocus);
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
        let tr = trace_full_ray(&din[i], &lens, refocus);
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
    let mut a: f64;
    let mut r: f64;
    let maxangle: f64 = 2.0 * PI;
    let mut count: usize = 0;

    let mut rng = rand::thread_rng();

    for _ in 0..gr.baserays
    {
        r = gr.half_ap * rng.gen_range(0.0_f64, 1.0).sqrt();
        a = rng.gen_range(0.0_f64, maxangle );
        x = r * a.cos();
        y = r * a.sin();

        let pvbase = Vector3D { x: x, y: y, z: 0.0 };

        for _ in 0..gr.num_angles
        {
            r = gr.half_ang  * rng.gen_range(0.0_f64, 1.0).sqrt();
            a = rng.gen_range(0.0_f64, maxangle );
            xdir = r * a.cos();
            ydir = r * a.sin();

            let edir = Vector3D {
                x: xdir,
                y: ydir,
                z: (1.0 - xdir * xdir - ydir * ydir).sqrt(),
            };
            din[count].pvector = pvbase.clone();
            din[count].edir = edir;
            count += 1;
        }
    }
}

pub fn gen_random_rays_while(gr: GenRays, din: &mut [Ray])
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

        for _ in 0..gr.num_angles
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
            din[count].pvector = pvbase.clone();
            din[count].edir = edir;
            count += 1;
        }
    }
}

#[no_mangle]
pub extern "C" fn trace_ray(p0: &Vector3D, e0: &Vector3D, lens: &Lens, refocus: f64) -> Ray
{
    // Trace ray from srf 0 to first lens surface. The axial distance here should be zero.
    let p2 = translate_to_surface(p0, e0, &lens.side1, 0.0);
    let n2 = calc_slope(&p2, &lens.side1);
    let e2 = calc_dir_sines(&e0, &n2, 1.0, lens.n); // after refraction

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
    
    //let one: f64 = p0.y;
    //return Ray{pvector: Vector3D{x: one, y: one, z: one}, edir: Vector3D{x: one, y: one, z: one}};
}


#[cfg(test)]
mod tests
{
    use super::*;

    #[test]
    fn it_works()
    {
        assert_eq!(2, 2)
    }
}
