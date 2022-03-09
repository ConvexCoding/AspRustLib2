
// ***************** FFT Tools ************************** */

use num_complex::Complex;


pub fn multi_conjugate_return_real(q: Vec<Vec<Complex<f64>>>) -> Vec<Vec<f64>> {

    let rows = q.len();
    let cols = q[0].len();

    let mut t = gen_zero_2d(rows);  
    for r in 0..rows {
        for c in 0..cols { 
            t[r][c] = (q[r][c] * q[r][c].conj()).re;
        }
    }
    return t;
}


pub fn _multi_conjugate(q: Vec<Vec<Complex<f64>>>) -> Vec<Vec<Complex<f64>>> {

    let rows = q.len();
    let cols = q[0].len();

    let mut t = gen_zero_complex_2d(rows);  
    for r in 0..rows {
        for c in 0..cols { 
            t[r][c] = q[r][c] * q[r][c].conj();
        }
    }
    return t;
}


pub fn _find_min_max(data: &Vec<Vec<f64>>) -> (f64, f64) {

    let rows = data.len();
    let cols = data[0].len();

    let mut min: f64 = 1.0e20;
    let mut max: f64 = -1.0e20;

    for r in 0..rows {
        for c in 0..cols { 
            if data[r][c] > max {
                max = data[r][c];
            }
            if data[r][c] < min {
                min = data[r][c];
            }
        }
    }
    return (min, max);

}


pub fn find_max(data: Vec<Vec<f64>>) -> f64 {

    let rows = data.len();
    let cols = data[0].len();

    let mut max: f64 = -1.0e20;

    for r in 0..rows {
        for c in 0..cols { 
            if data[r][c] > max {
                max = data[r][c];
            }
        }
    }
    return max;
}

pub fn _find_max_complex(data: Vec<Vec<Complex<f64>>>) -> f64 {

    let rows = data.len();
    let cols = data[0].len();

    let mut max: f64 = -1.0e20;

    for r in 0..rows {
        for c in 0..cols { 
            if data[r][c].re > max {
                max = data[r][c].re;
            }
        }
    }
    return max;
}

pub fn scaledata(data: &mut Vec<Vec<f64>>, xfactor: f64) {

    let rows = data.len();
    let cols = data[0].len();

    for r in 0..rows {
        for c in 0..cols { 
            data[r][c] /= xfactor;
        }
    }
}


pub fn gen_zero_complex_2d(numberelements: usize) -> Vec<Vec<Complex<f64>>> {
  
    let mut data: Vec<Vec<Complex<f64>>> = vec!();
    let row: Vec<Complex<f64>> = vec!(Complex{re: 0.0, im: 0.0}; numberelements);
    for _ in 0..numberelements {
        data.push(row.clone());
    }
    return data;
}


pub fn gen_zero_2d(numberelements: usize) -> Vec<Vec<f64>> {
    let mut data: Vec<Vec<f64>> = vec!();
    let row: Vec<f64> = vec!(0.0; numberelements);
    for _ in 0..numberelements {
        data.push(row.clone());
    }
    return data;
}


pub fn intlog2(x: u32) -> u32 {
    return (x as f64).log2().ceil() as u32;
}


pub fn _ispowerof2( x: u32 ) -> bool {
    if ( x > 0 ) && ( ( x & ( x - 1 ) ) == 0 ) {
        return true;
    } else {
        return false;
    }
}


pub fn _power2(x: u32) -> u32 {
    if x > 0 && x < 31 {
        return 1 << x;
    }    
    else {
        return 0;
    }
}

pub fn get_complex_vec(wfe: &Vec<Vec<f64>>, mask: &Vec<Vec<f64>>, padsize: usize ) -> Vec<Vec<Complex<f64>>> {

    let ix = mask.len();
    let iy = mask[0].len();

    let mut data: Vec<Vec<Complex<f64>>> = gen_zero_complex_2d(padsize);

    let start = (padsize - mask[0].len()) / 2;

    for r in 0..ix {
        for c in 0..iy {
            data[start + r][start + c] = mask[r][c] * Complex::exp(Complex{re: 0.0, im: wfe[r][c]});
        }
    }

    return data;
}

pub fn slicecore(data: Vec<Vec<f64>>, core: usize) -> Vec<Vec<f64>> {
    let rows = data.len();
    let cols = data[0].len();

    let startrow = (rows - core) / 2;
    let startcol = (cols - core) / 2;

    let mut odata = gen_zero_2d(core);

    for r in 0..core {
        for c in 0..core {
            odata[r][c] = data[startrow + r][startcol + c];
        }
    }

    return odata;
}