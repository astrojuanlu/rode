import numpy as np
import plotly.graph_objects as go

from rode import plate_displacement_field_py


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

    ww = np.zeros_like(xx)

    plate_displacement_field_py(
        xx,
        yy,
        ww,
        xi,
        eta,
        P,
        D,
        l_x,
        l_y,
        max_m,
        max_n
    )

    fig = go.Figure()
    fig.add_surface(z=ww)
    fig.show()
