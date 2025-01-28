from math import pi

from rode import euler_method_demo_py


def y(x):
    return -x


if __name__ == "__main__":
    y0_demo = 1.0
    t_start_demo = 0.0
    t_end_demo = 2.0 * pi
    num_points_demo = 100
    results = euler_method_demo_py(y, y0_demo, t_start_demo, t_end_demo, num_points_demo)

    print(type(results))
    print(results)
