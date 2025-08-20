#!/usr/bin/env python3
"""
Quick verification test for all three plate displacement implementations.
"""

import numpy as np
import time

from rode import plate_displacement_field_py, plate_displacement_field_par_py
from rode.navier import plate_displacement_field_py as plate_displacement_field_py_pure


def test_implementations():
    # Test with different parameters to ensure robustness
    l_x, l_y = 2.0, 1.5  # Non-square plate
    h = 25e-3
    E = 200e9  # Steel
    nu = 0.3
    max_m, max_n = 12, 10  # Different series terms
    NUM_POINTS = 51  # Smaller grid for faster testing
    P = -5e3
    xi, eta = 0.7, 0.6  # Off-center load

    D = h**3 * E / (12 * (1 - nu**2))

    x = np.linspace(0, l_x, num=NUM_POINTS)
    y = np.linspace(0, l_y, num=NUM_POINTS)
    xx, yy = np.meshgrid(x, y)

    print(f"Testing with {NUM_POINTS}x{NUM_POINTS} grid, max_m={max_m}, max_n={max_n}")
    print(f"Plate dimensions: {l_x}x{l_y}m, Load at ({xi}, {eta})")
    print()

    # Test all three versions
    ww_python = np.zeros_like(xx)
    start = time.time()
    plate_displacement_field_py_pure(
        xx, yy, ww_python, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    t_python = time.time() - start

    ww_rust_serial = np.zeros_like(xx)
    start = time.time()
    plate_displacement_field_py(
        xx, yy, ww_rust_serial, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    t_rust_serial = time.time() - start

    ww_rust_parallel = np.zeros_like(xx)
    start = time.time()
    plate_displacement_field_par_py(
        xx, yy, ww_rust_parallel, xi, eta, P, D, l_x, l_y, max_m, max_n
    )
    t_rust_parallel = time.time() - start

    # Check differences
    diff_py_rs = np.max(np.abs(ww_python - ww_rust_serial))
    diff_rs_rp = np.max(np.abs(ww_rust_serial - ww_rust_parallel))

    print(
        f"Times: Python={t_python:.4f}s, Rust serial={t_rust_serial:.4f}s, Rust parallel={t_rust_parallel:.4f}s"
    )
    print(f"Max differences: Py-Rs={diff_py_rs:.2e}, Rs-Rp={diff_rs_rp:.2e}")
    print(
        f"Speedups: Rust serial {t_python / t_rust_serial:.1f}x, Rust parallel {t_python / t_rust_parallel:.1f}x"
    )
    print(f"Max displacement: {np.max(np.abs(ww_rust_serial)):.3e} m")

    # Verify correctness
    if diff_py_rs < 1e-13 and diff_rs_rp < 1e-15:
        print("✅ All implementations produce identical results!")
        return True
    else:
        print("❌ Results differ!")
        return False


if __name__ == "__main__":
    success = test_implementations()
    exit(0 if success else 1)
