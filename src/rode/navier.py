import numpy as np
from numpy import pi, sin

from numba import njit as jit


@jit
def plate_displacement_point(
    x: float,
    y: float,
    xi: float,
    eta: float,
    p_load: float,
    d_rigidity: float,
    l_x: float,
    l_y: float,
    m_max: int,
    n_max: int,
) -> float:
    """Compute plate displacement at a single point using Navier series.

    This implementation matches the Rust version exactly.
    """
    result = 0.0

    for mm in range(1, m_max + 1):
        for nn in range(1, n_max + 1):
            result += (
                sin(mm * pi * x / l_x)
                * sin(nn * pi * y / l_y)
                * sin(mm * pi * xi / l_x)
                * sin(nn * pi * eta / l_y)
            ) / ((mm / l_x) ** 2 + (nn / l_y) ** 2) ** 2

    return result * 4.0 * p_load / (pi**4 * d_rigidity * l_x * l_y)


@jit
def plate_displacement_field_py(
    xx: np.ndarray,
    yy: np.ndarray,
    ww: np.ndarray,
    xi: float,
    eta: float,
    p_load: float,
    d_rigidity: float,
    l_x: float,
    l_y: float,
    m_max: int,
    n_max: int,
) -> None:
    """Compute plate displacement field using Navier series.

    This implementation matches the Rust version exactly.
    """
    max_i, max_j = ww.shape

    for ii in range(max_i):
        for jj in range(max_j):
            ww[ii, jj] = plate_displacement_point(
                xx[ii, jj],
                yy[ii, jj],
                xi,
                eta,
                p_load,
                d_rigidity,
                l_x,
                l_y,
                m_max,
                n_max,
            )


# Legacy function for backwards compatibility (kept as-is)
@jit
def a_mn_point(P, a, b, xi, eta, mm, nn):
    """Navier series coefficient for concentrated load."""
    return 4 * P * sin(mm * pi * xi / a) * sin(nn * pi * eta / b) / (a * b)


@jit
def plate_displacement(xx, yy, ww, a, b, P, xi, eta, D, max_m, max_n):
    max_i, max_j = ww.shape
    for mm in range(1, max_m):
        for nn in range(1, max_n):
            for ii in range(max_i):
                for jj in range(max_j):
                    a_mn = a_mn_point(P, a, b, xi, eta, mm, nn)
                    ww[ii, jj] += (
                        a_mn
                        / (mm**2 / a**2 + nn**2 / b**2) ** 2
                        * sin(mm * pi * xx[ii, jj] / a)
                        * sin(nn * pi * yy[ii, jj] / b)
                        / (pi**4 * D)
                    )
