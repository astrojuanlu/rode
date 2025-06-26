import numpy as np
import plotly.graph_objects as go

from plate import plate_displacement_field


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

    plate_displacement_field(
        xx,
        yy,
        ww,
        P,
        xi,
        eta,
        l_x,
        l_y,
        h,
        E,
        nu,
        D,
        max_m=max_m,
        max_n=max_n
    )

    fig = go.Figure()
    fig.add_surface(z=ww)
    fig.show()
