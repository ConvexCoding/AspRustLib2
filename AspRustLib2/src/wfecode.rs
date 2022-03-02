use crate::ray_vector::WFE_Ray;

#[no_mangle]
extern "C" fn tracewferays(
    ptr: *mut WFE_Ray,
    npts: usize,
    loopsize: i32,
    lens: Lens,
    refocus: f64,
) -> i32
{
    let din: &mut [WFE_Ray] = unsafe { std::slice::from_raw_parts_mut(ptr, npts) };

    gen_wfe_rays(5.0f64, npts, din, loopsize);

    npts as i32
}

fn gen_wfe_rays(apert: f64, size: usize, din: &mut [WFE_Ray], loopsize: i32)
{
    let mut x: f64;
    let mut y: f64;
    let diag = apert * apert * 1.00001f64;

    let step = 2f64 * apert / (loopsize - 1i32) as f64;
    let mut count = 0usize;

    let dirvector = Vector3D {
        x: 0.0f64,
        y: 0.0f64,
        z: 1.0f64,
    };

    for row in 0..loopsize
    {
        y = apert - row as f64 * step;
        for col in 0..loopsize
        {
            x = -apert + col as f64 * step;

            din[count].rstart.pvector = Vector3D { x: x, y: y, z: 0.0 };
            din[count].rstart.edir = dirvector;
            din[count].ix = col;
            din[count].iy = row;
            din[count].isvalid = diag > x * x + y * y;

            count = count + 1i32;
        }
    }
}

fn calc_wfe(wferay: WFE_Ray, lens: &Lens, refocus: f64)
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
    let rsqsq = rsq * rsq;

    let rm = trace_ray(p0, e0, lens, 0.0f64);
    //let ym = rm.pvector;
    let (ymaoi, ymlsa) = calcstuff(rm);

    let rz = trace_ray(p1, e0, lens, 0.0f64);
    //let yz = rz.pvector;
    let (_yzaoi, yzlsa) = calcstuff(rz);

    wferay.rend = trace_ray(p0, e0, lens, refocus);
    //let yfinal = rfinal.pvector;
    let (yfaoi, yflsa) = calcstuff(rfinal);
    wferay.lsa = yflsa;

    let a = (4.0f64 * yzlsa - ymlsa) / rsq;
    let b = (2.0f64 * ymlsa - 4.0f64 * yzlsa) / rsqsq;

    wferay.opd = 1000.0f64
        * (ymaoi.sin() * ymaoi.sin() / 2.0f64)
        * (refocus - a * rsq / 2.0f64 - b * rsqsq / 3.0f64)
        / lens.wl;
}

pub fn calcstuff(ray: Ray) -> (f64, f64)
{
    let eflat = Vector3D {
        x: 0.0f64,
        y: 0.0f64,
        z: 1.0f64,
    };
    let edot = dot_product(ray.edir, eflat);
    let aoi = edot.acos();
    let lsa = -1.0f64 * (ray.pvector.x * ray.pvector.x + ray.pvector.y * ray.pvector.y).sqrt()
        / aoi.tan();
    //let lsa = -ray.pvector.ydir / ray.edir.y;
    (aoi, lsa)
}

fn calc_wfe(wferay: &mut WFE_Ray, lens: &Lens, refocus: f64)
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
    let rsqsq = rsq * rsq;

    let rm = trace_ray(p0, e0, lens, 0.0f64);
    //let ym = rm.pvector;
    let (ymaoi, ymlsa) = calcstuff(rm);

    let rz = trace_ray(p1, e0, lens, 0.0f64);
    //let yz = rz.pvector;
    let (_yzaoi, yzlsa) = calcstuff(rz);

    let rfinal = trace_ray(p0, e0, lens, refocus);
    //let yfinal = rfinal.pvector;
    let (yfaoi, yflsa) = calcstuff(rfinal);
    wferay.lsa = yflsa;

    let a = (4.0f64 * yzlsa - ymlsa) / rsq;
    let b = (2.0f64 * ymlsa - 4.0f64 * yzlsa) / rsqsq;

    wferay.opd = 1000.0f64
        * (ymaoi.sin() * ymaoi.sin() / 2.0f64)
        * (refocus - a * rsq / 2.0f64 - b * rsqsq / 3.0f64)
        / lens.wl;
}
