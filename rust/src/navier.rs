use std::f64::consts::PI;

use numpy::ndarray::{Array2, ArrayViewMut2};
use rayon::prelude::*;

pub fn plate_displacement(
    x: f64,
    y: f64,
    xi: f64,
    eta: f64,
    p_load: f64,
    d_rigidity: f64,
    l_x: f64,
    l_y: f64,
    m_max: i32,
    n_max: i32,
) -> f64 {
    let mut result = 0.0;

    for mm in 1..=m_max {
        for nn in 1..=n_max {
            result += ((mm as f64 * PI * x / l_x).sin()
                * (nn as f64 * PI * y / l_y).sin()
                * (mm as f64 * PI * xi / l_x).sin()
                * (nn as f64 * PI * eta / l_y).sin())
                / ((mm as f64 / l_x).powi(2) + (nn as f64 / l_y).powi(2)).powi(2)
        }
    }
    result * 4.0 * p_load / (PI.powi(4) * d_rigidity * l_x * l_y)
}

pub fn plate_displacement_field(
    xx: &Array2<f64>,
    yy: &Array2<f64>,
    ww: &mut ArrayViewMut2<f64>,
    xi: f64,
    eta: f64,
    p_load: f64,
    d_rigidity: f64,
    l_x: f64,
    l_y: f64,
    m_max: i32,
    n_max: i32,
) -> () {
    let shape = ww.dim();
    let max_i = shape.0;
    let max_j = shape.1;

    for ii in 0..max_i {
        for jj in 0..max_j {
            ww[(ii, jj)] = plate_displacement(
                xx[(ii, jj)],
                yy[(ii, jj)],
                xi,
                eta,
                p_load,
                d_rigidity,
                l_x,
                l_y,
                m_max,
                n_max,
            );
        }
    }
}

pub fn plate_displacement_field_par(
    xx: &Array2<f64>,
    yy: &Array2<f64>,
    ww: &mut ArrayViewMut2<f64>,
    xi: f64,
    eta: f64,
    p_load: f64,
    d_rigidity: f64,
    l_x: f64,
    l_y: f64,
    m_max: i32,
    n_max: i32,
) -> () {
    let shape = ww.dim();
    let max_i = shape.0;
    let max_j = shape.1;

    // Create index pairs for parallel iteration
    let indices: Vec<(usize, usize)> = (0..max_i)
        .flat_map(|ii| (0..max_j).map(move |jj| (ii, jj)))
        .collect();

    // Collect results from parallel computation
    let results: Vec<f64> = indices
        .par_iter()
        .map(|(ii, jj)| {
            plate_displacement(
                xx[(*ii, *jj)],
                yy[(*ii, *jj)],
                xi,
                eta,
                p_load,
                d_rigidity,
                l_x,
                l_y,
                m_max,
                n_max,
            )
        })
        .collect();

    // Write results back to the output array
    for (idx, (ii, jj)) in indices.iter().enumerate() {
        ww[(*ii, *jj)] = results[idx];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use numpy::ndarray::Array;

    #[test]
    fn test_plate_displacement_field_par_matches_serial() {
        // Test parameters based on plate.py
        let l_x = 1.0;
        let l_y = 1.0;
        let h = 50e-3_f64;
        let e = 69e9_f64;
        let nu = 0.35_f64;
        let max_m = 8; // Smaller for testing
        let max_n = 8;
        let p = -10e3_f64;
        let xi = l_x / 2.0;
        let eta = l_y / 2.0;
        let d = h.powi(3) * e / (12.0 * (1.0 - nu.powi(2)));

        // Create a small test grid
        let num_points = 11;
        let x_vals: Vec<f64> = (0..num_points)
            .map(|i| (i as f64) * l_x / ((num_points - 1) as f64))
            .collect();
        let y_vals: Vec<f64> = (0..num_points)
            .map(|i| (i as f64) * l_y / ((num_points - 1) as f64))
            .collect();

        // Create meshgrid
        let mut xx_data = Vec::new();
        let mut yy_data = Vec::new();
        for y in &y_vals {
            for x in &x_vals {
                xx_data.push(*x);
                yy_data.push(*y);
            }
        }

        let xx = Array::from_shape_vec((num_points, num_points), xx_data).unwrap();
        let yy = Array::from_shape_vec((num_points, num_points), yy_data).unwrap();

        // Create output arrays
        let mut ww_serial = Array::zeros((num_points, num_points));
        let mut ww_parallel = Array::zeros((num_points, num_points));

        // Test serial version
        plate_displacement_field(
            &xx,
            &yy,
            &mut ww_serial.view_mut(),
            xi,
            eta,
            p,
            d,
            l_x,
            l_y,
            max_m,
            max_n,
        );

        // Test parallel version
        plate_displacement_field_par(
            &xx,
            &yy,
            &mut ww_parallel.view_mut(),
            xi,
            eta,
            p,
            d,
            l_x,
            l_y,
            max_m,
            max_n,
        );

        // Compare results - they should be identical
        for ii in 0..num_points {
            for jj in 0..num_points {
                let diff = (ww_serial[(ii, jj)] - ww_parallel[(ii, jj)]).abs();
                assert!(
                    diff < 1e-15,
                    "Results differ at ({}, {}): serial={}, parallel={}, diff={}",
                    ii,
                    jj,
                    ww_serial[(ii, jj)],
                    ww_parallel[(ii, jj)],
                    diff
                );
            }
        }

        // Also check that we get non-zero results
        let max_displacement = ww_serial
            .iter()
            .map(|x| x.abs())
            .fold(0.0_f64, |a, b| a.max(b));
        assert!(max_displacement > 0.0, "Displacement should be non-zero");
    }

    #[test]
    fn test_plate_displacement_single_point() {
        // Test a single point calculation
        let x = 0.5;
        let y = 0.5;
        let xi = 0.5;
        let eta = 0.5;
        let p_load = -1000.0;
        let d_rigidity = 1000.0;
        let l_x = 1.0;
        let l_y = 1.0;
        let m_max = 10;
        let n_max = 10;

        let displacement = plate_displacement(
            x, y, xi, eta, p_load, d_rigidity, l_x, l_y, m_max, n_max,
        );

        // The displacement should be negative for a negative load
        assert!(displacement < 0.0, "Displacement should be negative for negative load");
    }
}
