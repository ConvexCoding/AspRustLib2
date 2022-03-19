use std::ptr;
use super::lens_struct::{Lens, OptVariables};
use super::ray_vector::{Ray, Vector3D, CPROPV};
use super::rtm_subs::trace_full_ray;

//pub fn opti_lens(rays: &[f64], lensin: Lens, optv: &OptVariables) -> Lens
pub fn opti_lens(ray_ys: &[f64], lensin: Lens, optv: &OptVariables) -> Lens
{
    let mut lens = lensin;

    let mut rays: Vec<Ray>  = Vec::with_capacity(ray_ys.len());
    for i in 0..ray_ys.len(){
        rays.push( Ray{pvector: Vector3D{x: 0.0, y: ray_ys[i] * lens.ap, z: 0.0}, edir: CPROPV});
    }
    
    let (mut cc, mut ad, mut ae);
    let loopsize = 15;

    for _k in 0..loopsize 
    {
        for _j in 0..loopsize
        {
            let (mut deltk, mut deltad, mut deltae, mintweakstep) = (0.6, 3.1e-7, 1.1e-9, 1e-15);

            for _i in 0..loopsize
            {                 
                if optv.iscc1on == 1
                {
                    cc = finddirection(&mut lens, &rays, deltk, 0);
                    lens.side1.k += cc;
                    if cc.abs() < mintweakstep {
                        deltk /= 10.0;
                    }
                }

                if optv.isad1on == 1
                {
                    ad = finddirection(&mut lens, &rays, deltad, 1);
                    lens.side1.ad += ad;
                    if ad.abs() < mintweakstep {
                        deltad /= 10.0;
                    }
                }

                if optv.isae1on == 1
                {
                    ae = finddirection(&mut lens, &rays, deltae, 2);
                    lens.side1.ae += ae;
                    if ae.abs() < mintweakstep {
                        deltae /= 10.0;
                    }
                }

                if optv.iscc2on == 1
                {
                    cc = finddirection(&mut lens, &rays, deltk, 3);
                    lens.side2.k += cc;
                    if cc.abs() < mintweakstep {
                        deltk /= 10.0;
                    }
                }

                if optv.isad2on == 1
                {
                    ad = finddirection(&mut lens, &rays, deltad, 4);
                    lens.side2.ad += ad;
                    if ad.abs() < mintweakstep {
                        deltad /= 10.0;
                    }
                }

                if optv.isae2on == 1
                {
                    ae = finddirection(&mut lens, &rays, deltae, 5);
                    lens.side2.ae += ae;
                    if ae.abs() < mintweakstep {
                        deltae /= 10.0;
                    }
                }
            }
        }
    }

    return lens;
}


pub fn finddirection(lens: &mut Lens, rays: &Vec<Ray>, delta_x: f64, whichvar: i32) -> f64
{
    let raw_ptr: *mut f64 = match whichvar
    {
        0 => (&mut lens.side1.k),
        1 => (&mut lens.side1.ad),
        2 => (&mut lens.side1.ae),
        3 => (&mut lens.side2.k),
        4 => (&mut lens.side2.ad),
        5 => (&mut lens.side2.ae),
        _ => (ptr::null_mut())
    };
    
    if raw_ptr.is_null() {
        return 0.0;
    };

    /* 
    let raw_ptr: *mut f64;
    if whichvar == 0 {
        raw_ptr = &mut lens.side1.k;
    } else if whichvar == 1 {
        raw_ptr = &mut lens.side1.ad;
    } else if whichvar == 2 {
        raw_ptr = &mut lens.side1.ae;
    } else if whichvar == 3 {
        raw_ptr = &mut lens.side2.k;
    } else if whichvar == 4 {
        raw_ptr = &mut lens.side2.ad;
    } else if whichvar == 5 {
        raw_ptr = &mut lens.side2.ae;
    } else {
        return 0.0;
    }
    */

    let ( basevalue, center, right, left);

    unsafe { 
        basevalue = *raw_ptr;
        center = calc_err_slim(&lens, &rays);
        *raw_ptr += delta_x;
        right = calc_err_slim(&lens, &rays);
        *raw_ptr -= 2.0_f64 * delta_x;
        left = calc_err_slim(&lens, &rays);
        *raw_ptr = basevalue;
    }

    if left < center { 
        return -delta_x;
    }
    if right < center {
        return delta_x;
    }
    return 0.0;
}

pub fn calc_err_slim(lens: &Lens, rays: &Vec<Ray>) -> f64
{
    let rayct = rays.len();
    let mut sumsum: f64 = 0.0;

    for i in 0..rayct{
        let rayout = trace_full_ray(&rays[i], &lens, 0.0_f64);
        sumsum += rayout.pvector.y * rayout.pvector.y;
    }
    return (sumsum / rayct as f64).sqrt();
}

