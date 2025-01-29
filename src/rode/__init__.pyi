import typing as t
import optype.numpy as onpt
import numpy as np

def euler_method_demo_py(
    dx_dt: t.Callable[[float], float],
    y0: float,
    t_start: float,
    t_end: float,
    num_points: int,
) -> onpt.Array1D[np.float]: ...

def euler_method_demo_alt_py(
    dx_dt: t.Callable[[float], float],
    y0: float,
    t: np.ndarray,
    h: float,
) -> onpt.Array1D[np.float]: ...

def farnocchia_coe_py(
    k_kms: float,
    p_km: float,
    ecc: float,
    inc_rad: float,
    raan_rad: float,
    argp_rad: float,
    nu_rad: float,
    tof_s: float,
) -> float: ...
