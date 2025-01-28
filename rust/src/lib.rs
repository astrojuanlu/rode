use numpy::ndarray::Array1;
use numpy::{IntoPyArray, PyArray1};
use pyo3::prelude::*;

/// Implements the Euler method to solve ordinary differential equations.
pub fn euler_method<F>(dx_dt: F, y0: f64, t: &Array1<f64>, h: f64) -> Array1<f64>
where
    F: Fn(f64) -> f64,
{
    let mut y = Array1::zeros(t.len());
    y[0] = y0;
    for i in 1..t.len() {
        y[i] = y[i - 1] + dx_dt(y[i - 1]) * h;
    }
    y
}

/// A demo function that calls euler_method with a dummy function.
pub fn euler_method_demo<F>(
    dx_dt: F,
    y0: f64,
    t_start: f64,
    t_end: f64,
    num_points: usize,
) -> Array1<f64>
where
    F: Fn(f64) -> f64,
{
    let t = Array1::linspace(t_start, t_end, num_points);
    let h = t[1] - t[0];
    euler_method(dx_dt, y0, &t, h)
}

/// A wrapper function to expose euler_method_demo to Python.
#[pyfunction]
fn euler_method_demo_py(
    py: Python,
    y0: f64,
    t_start: f64,
    t_end: f64,
    num_points: usize,
) -> PyResult<Bound<'_, PyArray1<f64>>> {
    let dx_dt = |x: f64| -> f64 { -x };
    let result = euler_method_demo(dx_dt, y0, t_start, t_end, num_points);
    Ok(result.into_pyarray(py))
}

/// A Python module implemented in Rust.
#[pymodule]
#[pyo3(name = "_rode")]
fn rode(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(euler_method_demo_py, m)?)?;
    Ok(())
}
