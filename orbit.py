import numpy as np

from rode import farnocchia_coe_py


if __name__ == "__main__":
    k = 398600.4418
    p = 6780.84
    ecc = 0.00130547
    inc = np.deg2rad(51.6012)
    raan = np.deg2rad(198.3795)
    argp = np.deg2rad(39.2629)
    nu = np.deg2rad(46.5958)

    tof = 20_000.0

    nu = farnocchia_coe_py(k, p, ecc, inc, raan, argp, nu, tof)
    print(f"nu = {np.rad2deg(nu)} deg")
