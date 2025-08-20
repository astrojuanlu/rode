import numpy as np
import plotly.graph_objects as go
import time

from rode import plate_displacement_field_py, plate_displacement_field_par_py


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

    # Test serial version
    ww_serial = np.zeros_like(xx)

    start_time = time.time()
    plate_displacement_field_py(
        xx, yy, ww_serial, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    serial_time = time.time() - start_time

    # Test parallel version
    ww_parallel = np.zeros_like(xx)

    start_time = time.time()
    plate_displacement_field_par_py(
        xx, yy, ww_parallel, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    parallel_time = time.time() - start_time

    # Verify results are identical
    max_diff = np.max(np.abs(ww_serial - ww_parallel))
    print(f"Maximum difference between serial and parallel: {max_diff:.2e}")

    if max_diff < 1e-15:
        print("✅ Results are identical!")
    else:
        print("❌ Results differ!")

    # Performance comparison
    print(f"Serial time: {serial_time:.4f} seconds")
    print(f"Parallel time: {parallel_time:.4f} seconds")
    if parallel_time > 0:
        speedup = serial_time / parallel_time
        print(f"Speedup: {speedup:.2f}x")

    # Maximum displacement
    max_w = np.max(np.abs(ww_serial))
    print(f"Maximum displacement: {max_w:.3e} m")

    fig = go.Figure()
    fig.add_surface(z=ww_serial, name="Serial")
    fig.update_layout(
        title="Plate Displacement (Serial and Parallel versions are identical)"
    )
    fig.show(renderer="browser")
