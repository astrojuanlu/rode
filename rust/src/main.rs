use numpy::ndarray::Array1;
use std::f64::consts::PI;

use rode::euler_method;

fn main() {
    // Demonstrate euler_method
    let dx_dt = |x: f64| -> f64 { -x };
    let y0 = 1.0;
    let t = Array1::linspace(0.0, 2.0 * PI, 100);
    let h = t[1] - t[0];

    let result = euler_method(dx_dt, y0, &t, h);

    println!("Results from euler_method:");
    for (i, &value) in result.iter().enumerate() {
        println!("t = {:.2}, y = {:.5}", t[i], value);
    }
}
