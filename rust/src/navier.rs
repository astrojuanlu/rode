use std::f64::consts::PI;

use numpy::ndarray::{Array2, ArrayViewMut2};

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
