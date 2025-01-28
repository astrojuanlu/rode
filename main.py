import numpy as np

from rode import euler_method_demo_py, euler_method_demo_alt_py


def y(x):
    return -x


if __name__ == "__main__":
    y0_demo = 1.0
    t_start = 0.0
    t_end = 2.0 * np.pi
    num_points = 100
    results = euler_method_demo_py(y, y0_demo, t_start, t_end, num_points)

    print(type(results))
    print(results)

    t = np.linspace(t_start, t_end, num_points)
    h = t[1] - t[0]

    results = euler_method_demo_alt_py(y, y0_demo, t, h)

    print(type(results))
    print(results)
