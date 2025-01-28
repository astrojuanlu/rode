use numpy::ndarray::Array1;

/// Euler method to solve ordinary differential equations
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

/// Alternative Euler method with time bounds and fixed number of points
pub fn euler_method_alt<F>(
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
