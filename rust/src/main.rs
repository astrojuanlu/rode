use ndarray::Array1;
use std::f64::consts::PI;

use rode::{euler_method, euler_method_demo};

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

    // Demonstrate euler_method_demo
    let y0_demo = 1.0;
    let t_start_demo = 0.0;
    let t_end_demo = 2.0 * PI;
    let num_points_demo = 100;

    let result_demo = euler_method_demo(y0_demo, t_start_demo, t_end_demo, num_points_demo);

    println!("\nResults from euler_method_demo:");
    let t_demo = Array1::linspace(t_start_demo, t_end_demo, num_points_demo);
    for (i, &value) in result_demo.iter().enumerate() {
        println!("t = {:.2}, y = {:.5}", t_demo[i], value);
    }
}
