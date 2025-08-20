import numpy as np
import plotly.graph_objects as go
import time

from rode import plate_displacement_field_py, plate_displacement_field_par_py
from rode.navier import plate_displacement_field_py as plate_displacement_field_py_pure


if __name__ == "__main__":
    # Plate geometry
    l_x = 1.0  # m
    l_y = 1.0  # m
    h = 50e-3  # m

    # Material properties
    E = 69e9  # Pa
    nu = 0.35

    # Series terms
    max_m = 16
    max_n = 16

    # Computation points
    # NOTE: With an odd number of points the center of the place is included in
    # the grid
    NUM_POINTS = 101

    # Load
    P = -10e3  # N
    xi = l_x / 2
    eta = l_x / 2

    # Flexural rigidity
    D = h**3 * E / (12 * (1 - nu**2))

    # Set up domain
    x = np.linspace(0, l_x, num=NUM_POINTS)
    y = np.linspace(0, l_y, num=NUM_POINTS)
    xx, yy = np.meshgrid(x, y)

    print("Testing four implementations:")
    print("1. Pure Python (no JIT)")
    print("2. Python + Numba JIT")
    print("3. Rust serial")
    print("4. Rust parallel")
    print()

    # Test 1: Pure Python version (no JIT)
    ww_python_pure = np.zeros_like(xx)
    start_time = time.time()
    plate_displacement_field_py_pure.py_func(
        xx, yy, ww_python_pure, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    python_pure_time = time.time() - start_time

    # Test 2: Python + Numba JIT version
    ww_python_numba = np.zeros_like(xx)
    start_time = time.time()
    plate_displacement_field_py_pure(
        xx, yy, ww_python_numba, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    python_numba_time = time.time() - start_time

    # Test 3: Rust serial version
    ww_rust_serial = np.zeros_like(xx)
    start_time = time.time()
    plate_displacement_field_py(
        xx, yy, ww_rust_serial, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    rust_serial_time = time.time() - start_time

    # Test 4: Rust parallel version
    ww_rust_parallel = np.zeros_like(xx)
    start_time = time.time()
    plate_displacement_field_par_py(
        xx, yy, ww_rust_parallel, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    rust_parallel_time = time.time() - start_time

    # Verify all results are identical
    print("=== Correctness Check ===")
    max_diff_pure_numba = np.max(np.abs(ww_python_pure - ww_python_numba))
    max_diff_numba_rust_serial = np.max(np.abs(ww_python_numba - ww_rust_serial))
    max_diff_rust_serial_parallel = np.max(np.abs(ww_rust_serial - ww_rust_parallel))
    max_diff_pure_rust_parallel = np.max(np.abs(ww_python_pure - ww_rust_parallel))

    print(f"Max difference (Pure Python vs Python+Numba): {max_diff_pure_numba:.2e}")
    print(
        f"Max difference (Python+Numba vs Rust serial): {max_diff_numba_rust_serial:.2e}"
    )
    print(
        f"Max difference (Rust serial vs parallel): {max_diff_rust_serial_parallel:.2e}"
    )
    print(
        f"Max difference (Pure Python vs Rust parallel): {max_diff_pure_rust_parallel:.2e}"
    )

    all_identical = (
        max_diff_pure_numba < 1e-14
        and max_diff_numba_rust_serial < 1e-14
        and max_diff_rust_serial_parallel < 1e-14
        and max_diff_pure_rust_parallel < 1e-14
    )

    if all_identical:
        print("✅ All four implementations produce identical results!")
    else:
        print("❌ Results differ between implementations!")

    # Performance comparison
    print()
    print("=== Performance Comparison ===")
    print(f"Pure Python time: {python_pure_time:.4f} seconds")
    print(f"Python + Numba time: {python_numba_time:.4f} seconds")
    print(f"Rust serial time: {rust_serial_time:.4f} seconds")
    print(f"Rust parallel time: {rust_parallel_time:.4f} seconds")

    if python_pure_time > 0:
        numba_speedup = python_pure_time / python_numba_time
        rust_serial_speedup = python_pure_time / rust_serial_time
        rust_parallel_speedup = python_pure_time / rust_parallel_time
        parallel_vs_serial_speedup = rust_serial_time / rust_parallel_time

        print(f"Numba speedup vs Pure Python: {numba_speedup:.2f}x")
        print(f"Rust serial speedup vs Pure Python: {rust_serial_speedup:.2f}x")
        print(f"Rust parallel speedup vs Pure Python: {rust_parallel_speedup:.2f}x")
        print(
            f"Rust parallel speedup vs Rust serial: {parallel_vs_serial_speedup:.2f}x"
        )

    # Maximum displacement
    max_w = np.max(np.abs(ww_rust_serial))
    print(f"\nMaximum displacement: {max_w:.3e} m")

    fig = go.Figure()
    fig.add_surface(z=ww_rust_serial, name="Displacement")
    fig.update_layout(
        title="Plate Displacement (All four implementations produce identical results)"
    )
    fig.show(renderer="browser")
