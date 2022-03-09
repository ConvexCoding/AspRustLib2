//****************** CREDITS & short Explaination***************************** */
// Source very roughly mimics AForge Math Library
// AForge.NET framework
// http://www.aforgenet.com/framework/
//
// Copyright © Andrew Kirillov, 2005-2009
// andrew.kirillov@aforgenet.com
//
// FFT idea from Exocortex.DSP library
// http://www.exocortex.org/dsp/
//
// get data is reading the wfe and mask grids, and then padding them to final size
// totalgrid size is required to be a power of 2:  64, 128, ... 1024, 2048, 4096.
// total grid sizes smaller than 128 produces coarse data.
// grid sizes >2048 seems to take an inordinate amount of time to calc in c#
// a rough rule of thumb is that total grid size should be 8 to 10 times larger
// than the initial wfe grid

use std::vec;
use std::f64::consts::PI;
//use std::fs::File;
//use std::io::{BufRead, BufReader};
//use std::io::Write;  
use num_complex::Complex;

use crate::tools::ffttools::{find_max, gen_zero_complex_2d, intlog2};
use crate::tools::ffttools::{scaledata, multi_conjugate_return_real};

const   _MIN_LENGTH: u32 = 2;
const   _MAX_LENGTH: u32 = 16384;
const   _MIN_BITS: u32 = 1;
const   _MAX_BITS: u32 = 14;

pub fn rustfft (data: &mut Vec<Vec<Complex<f64>>>, datadl: &mut Vec<Vec<Complex<f64>>> ) -> Vec<Vec<f64>>
{   
    let totalgrid: usize = data.len();
    let numbits = intlog2(totalgrid as u32);
    //println!("Rust FFT Total Grid: {},  Number of Bits: {}", totalgrid, numbits);

    // Run main fft2d routine
    fft2d(data, numbits, totalgrid);
    fft2d(datadl, numbits, totalgrid);

    let datashift = fft2shift(&data);
    let datadlshirt = fft2shift(&datadl);

    let mut dataout = multi_conjugate_return_real(datashift);
    let datadlaout = multi_conjugate_return_real(datadlshirt);

    let maxdl = find_max(datadlaout.clone());   
    //scaledata(&mut datadlaout, maxdl);
    scaledata(&mut dataout, maxdl);

    //let _maxdlnorm = find_max(datadlaout.clone());
    //let mut _max = find_max(dataout.clone());
    //_max = find_max(dataout.clone());


    return dataout;
}


fn fft2d(data: &mut Vec<Vec<Complex<f64>>>, numbits: u32, totalgrid: usize) {

    // step 1 - generate reverse bit settings for size of array
    let rbits: Vec<usize> = get_reversed_bits( numbits, totalgrid);

    // step 2 - generate a rotation matrix for data
    let mut rotation: Vec<Vec<Complex<f64>>> = vec!();    
    for i in 0..numbits {
        rotation.push(get_complex_rotation(i));
    }

    // step 3 - process data first rows then columns

    // process rows
    let mut row: Vec<Complex<f64>> = vec!(Complex{re: 0.0, im: 0.0}; totalgrid);

    for i in 0..totalgrid 
    {
        for j in 0..totalgrid {         // copy row
            row[j] = data[i][j];
        }
        fft(&mut row, &rbits, &rotation);          // transform it
        for j in 0..totalgrid {          //copy back
           data[i][j] = row[j];
        }
    }
    
    // process columns
    let mut col: Vec<Complex<f64>> = vec!(Complex{re: 0.0, im: 0.0}; totalgrid);

    for j in 0..totalgrid 
    {
        for i in 0..totalgrid {
            col[i] = data[i][j];
        }
        fft(&mut col, &rbits, &rotation);
        for i in 0..totalgrid {
            data[i][j] = col[i];
        }
    }
}

fn fft(data: &mut Vec<Complex<f64>>, rbits: &Vec<usize>, crotate: &Vec<Vec<Complex<f64>>>) {
    let n = data.len() as usize;
    let m = intlog2(n as u32) as usize;

    // reorder data first
    reorder_data(data, &rbits);

    // compute FFT
    let mut tn:usize = 1;
    let mut tm: usize;

    for k in 1..m+1
    {
        let rotation = &crotate[k - 1];

        tm = tn;
        tn <<= 1;

        for i in 0..tm
        {
            let t = rotation[i];

            //for (int even = i; even < n; even += tn)
            for even in (i..n).step_by(tn)
            {
                let odd = even + tm;
                let ce = data[even];
                let co = data[odd];

                let tr = co.re * t.re - co.im * t.im;
                let ti = co.re * t.im + co.im * t.re;

                data[even] += Complex{re: tr, im: ti};

                data[odd] = Complex{re: ce.re - tr, im: ce.im - ti};
            }
        }
    }

    for i in 0..n {
        data[i] /= n as f64;
    }
}

fn fft2shift(q: &Vec<Vec<Complex<f64>>>) -> Vec<Vec<Complex<f64>>> {
    let rows = q.len();
    let cols = q[0].len();

    let mut s = gen_zero_complex_2d(rows);      // shift columns
    let mut mid = cols / 2;
    for c in mid..cols 
    {
        for r in 0..rows
        {
            s[r][c - mid] = q[r][c];
            s[r][c] = q[r][c - mid];
        }
    }

    let mut t = gen_zero_complex_2d(rows);     // shift rows    
    mid = rows / 2;
    for r in mid..rows
    {
        for c in 0..rows
        {
            t[r - mid][c] = s[r][c];
            t[r][c] = s[r - mid][c];
        }
    }

    return t;
}

fn get_complex_rotation(numberofbits: u32) -> Vec<Complex<f64>> {
    let n = 1 << (numberofbits - 0);
    let mut ur = 1.0;
    let mut ui = 0.0;
    let angle = PI / (n as f64);
    let wr = angle.cos();
    let wi = angle.sin();

    let mut t: f64;

    let mut rotation: Vec<Complex<f64>> = vec!();
    for _ in 0..n 
    {
        rotation.push(Complex{re: ur, im: ui});
        t = ur * wi + ui * wr;
        ur = ur * wr - ui * wi;
        ui = t;
    }
    return rotation;
}

fn get_reversed_bits(numberofbits: u32, totalgridsize: usize ) -> Vec<usize> {

    // this should throw some type of error but not sure how to do it
    //if (numberofbits < minBits) || (numberofbits > maxBits) 
    //    return null;
    
    let mut rbits: Vec<usize> = vec!();

    // calculate the array
    for i in 0..totalgridsize
    {
        let mut oldbits = i;
        let mut newbits = 0;

        for _ in 0..numberofbits {
            newbits = newbits << 1 | oldbits & 1;
            oldbits = oldbits >> 1;
        }
        rbits.push(newbits);
    }
    return rbits;
}

pub fn reorder_data(data: &mut Vec<Complex<f64>>, rbits: &Vec<usize>) {
    let datalen = data.len();

    // check data length
    //if ( ( len < minLength ) || ( len > maxLength ) || ( !Tools.IsPowerOf2( len ) ) )
    //	throw new ArgumentException( "Incorrect data length." );


    for i in 0..datalen {
        let s = rbits[i];

        if s as usize > i {
            let t = data[i];
            data[i] = data[s];
            data[s] = t;
        }
    }
}



/*
// functions used for standalone version of ff2

fn read_data_from_file(fnamestr: &str) -> Vec<Vec<f64>> {
    let fnum = BufReader::new(File::open(fnamestr).unwrap());

    let mapdata: Vec<Vec<f64>> = fnum.lines()
            .map(|l| l.unwrap().split(", ")
            .map(|number| number.parse().unwrap())
            .collect())
            .collect();
    return mapdata;
}

fn get_diff_ltd_data(paddedsize: usize ) -> Vec<Vec<Complex<f64>>> {
    //C:/Users/glher/Desktop/fftdata/
    let maskdata = read_data_from_file("C:/Users/glher/Desktop/fftdata/mask.txt");

    let ix = maskdata.len();
    let iy = maskdata[0].len();

    let mut data: Vec<Vec<Complex<f64>>> = gen_zero_complex_2d(paddedsize);

    let start = (paddedsize - maskdata[0].len()) / 2;

    //println!("Mask Grid Size, x Amp Grid Size, Start Position: {0} x {1}, {2}", ix, iy, start);

    for r in 0..ix {
        for c in 0..iy {
            data[start + r][start + c] = maskdata[r][c] * Complex::exp(Complex{re: 0.0, im: maskdata[r][c]});
        }
    }

    return data;
}

fn get_data(paddedsize: usize ) -> Vec<Vec<Complex<f64>>> {
    //C:/Users/glher/Desktop/fftdata/
    let maskdata = read_data_from_file("C:/Users/glher/Desktop/fftdata/mask.txt");
    let wfedata = read_data_from_file("C:/Users/glher/Desktop/fftdata/wfe.txt");

    let ix = maskdata.len();
    let iy = maskdata[0].len();

    let mut data: Vec<Vec<Complex<f64>>> = gen_zero_complex_2d(paddedsize);

    let start = (paddedsize - maskdata[0].len()) / 2;

    println!("Mask Grid Size, x Amp Grid Size, Start Position: {0} x {1}, {2}", ix, iy, start);

    for r in 0..ix {
        for c in 0..iy {
            data[start + r][start + c] = maskdata[r][c] * Complex::exp(Complex{re: 0.0, im: wfedata[r][c]});
        }
    }

    return data;
}
*/


