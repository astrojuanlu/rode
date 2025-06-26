use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::prelude::*;
use pyo3::types::PyAny;

use euler::{euler_method, euler_method_alt};
use navier::{plate_displacement, plate_displacement_field};
use orbit::farnocchia_coe;

/// A wrapper function to expose euler_method_demo to Python.
///
/// # Arguments
///
/// * `dx_dt` - A function that computes the derivative of `y` with respect to `t`.
/// * `y0` - The initial value of `y`.
/// * `t` - An array of time points at which to compute the solution.
/// * `h` - The step size.
///
/// # Returns
///
/// An array of `y` values computed at each time point in `t`.
///
/// # Example
///
/// ```python
/// import numpy as np
/// from rode import euler_method_demo_py
///
/// def dx_dt(x):
///     return -2 * x
///
/// y0 = 1.0
/// t_start = 0.0
/// t_end = 5.0
/// num_points = 100
///
/// result = euler_method_demo_py(dx_dt, y0, t_start, t_end, num_points)
/// print(result)
/// ```
#[pyfunction]
fn euler_method_demo_py(
    py: Python,
    dx_dt: Py<PyAny>,
    y0: f64,
    t_start: f64,
    t_end: f64,
    num_points: usize,
) -> PyResult<Bound<'_, PyArray1<f64>>> {
    let dx_dt_closure = |x: f64| -> f64 {
        let args = (x,);
        dx_dt.call1(py, args).unwrap().extract(py).unwrap()
    };
    let result = euler_method_alt(dx_dt_closure, y0, t_start, t_end, num_points);
    Ok(result.into_pyarray(py))
}

/// Another wrapper function to expose euler_method_demo to Python.
#[pyfunction]
fn euler_method_demo_alt_py<'py>(
    py: Python<'py>,
    dx_dt: Py<PyAny>,
    y0: f64,
    t: PyReadonlyArray1<'py, f64>,
    h: f64,
) -> PyResult<Bound<'py, PyArray1<f64>>> {
    let dx_dt_closure = |x: f64| -> f64 {
        let args = (x,);
        dx_dt.call1(py, args).unwrap().extract(py).unwrap()
    };
    let t_array = t.as_array().to_owned();
    let result = euler_method(dx_dt_closure, y0, &t_array, h);
    Ok(result.into_pyarray(py))
}

#[pyfunction]
fn farnocchia_coe_py(
    k_kms: f64,
    p_km: f64,
    ecc: f64,
    inc_rad: f64,
    raan_rad: f64,
    argp_rad: f64,
    nu_rad: f64,
    tof_s: f64,
) -> PyResult<f64> {
    Ok(farnocchia_coe(
        k_kms, p_km, ecc, inc_rad, raan_rad, argp_rad, nu_rad, tof_s,
    ))
}

#[pyfunction]
fn plate_displacement_py(
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
) -> PyResult<f64> {
    Ok(plate_displacement(
        x, y, xi, eta, p_load, d_rigidity, l_x, l_y, m_max, n_max,
    ))
}

#[pyfunction]
fn plate_displacement_field_py<'py>(
    py: Python<'py>,
    xx: PyReadonlyArray2<'py, f64>,
    yy: PyReadonlyArray2<'py, f64>,
    ww: PyArray2<'py, f64>,
    xi: f64,
    eta: f64,
    p_load: f64,
    d_rigidity: f64,
    l_x: f64,
    l_y: f64,
    m_max: i32,
    n_max: i32,
) -> PyResult<()> {
    let xx_array = xx.as_array().to_owned();
    let yy_array = yy.as_array().to_owned();
    if xx_array.shape() != yy_array.shape() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Input arrays xx and yy must have the same shape.",
        ));
    }
    if ww.shape() != xx_array.shape() {
        return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(
            "Output array ww must have the same shape as xx and yy.",
        ));
    }
    let ww_array = ww.as_array_mut();

    plate_displacement_field(
        &xx_array,
        &yy_array,
        &mut ww_array,
        xi,
        eta,
        p_load,
        d_rigidity,
        l_x,
        l_y,
        m_max,
        n_max,
    );
    Ok(())
}

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "_rode")]
fn rode(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(euler_method_demo_py, m)?)?;
    m.add_function(wrap_pyfunction!(euler_method_demo_alt_py, m)?)?;
    m.add_function(wrap_pyfunction!(farnocchia_coe_py, m)?)?;
    m.add_function(wrap_pyfunction!(plate_displacement_py, m)?)?;
    m.add_function(wrap_pyfunction!(plate_displacement_field_py, m)?)?;
    Ok(())
}

pub mod euler;
pub mod navier;
pub mod orbit;
