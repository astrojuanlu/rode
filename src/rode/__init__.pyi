import typing as t
import optype.numpy as onpt
import numpy as np

def euler_method_demo_py(
    dx_dt: t.Callable[[float], float],
    y0: float,
    t_start: float,
    t_end: float,
    num_points: int,
) -> onpt.Array1D[np.float]:
    ...
