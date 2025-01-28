use std::f64::consts::PI;

fn kepler_equation(E: f64, M: f64, ecc: f64) -> f64 {
    E_to_M(E, ecc) - M
}

fn kepler_equation_prime(E: f64, _M: f64, ecc: f64) -> f64 {
    1.0 - ecc * E.cos()
}

fn kepler_equation_hyper(F: f64, M: f64, ecc: f64) -> f64 {
    F_to_M(F, ecc) - M
}

fn kepler_equation_prime_hyper(F: f64, _M: f64, ecc: f64) -> f64 {
    ecc * F.cosh() - 1.0
}

fn newton_elliptic(x0: f64, M: f64, ecc: f64, tol: f64, maxiter: usize) -> f64 {
    let mut p0 = x0;
    for _ in 0..maxiter {
        let fval = kepler_equation(p0, M, ecc);
        let fder = kepler_equation_prime(p0, M, ecc);
        let newton_step = fval / fder;
        let p = p0 - newton_step;
        if (p - p0).abs() < tol {
            return p;
        }
        p0 = p;
    }
    f64::NAN
}

fn newton_hyperbolic(x0: f64, M: f64, ecc: f64, tol: f64, maxiter: usize) -> f64 {
    let mut p0 = x0;
    for _ in 0..maxiter {
        let fval = kepler_equation_hyper(p0, M, ecc);
        let fder = kepler_equation_prime_hyper(p0, M, ecc);
        let newton_step = fval / fder;
        let p = p0 - newton_step;
        if (p - p0).abs() < tol {
            return p;
        }
        p0 = p;
    }
    f64::NAN
}

fn D_to_nu(D: f64) -> f64 {
    2.0 * D.atan()
}

fn nu_to_D(nu: f64) -> f64 {
    nu.tan() / 2.0
}

fn nu_to_E(nu: f64, ecc: f64) -> f64 {
    2.0 * ((1.0 - ecc) / (1.0 + ecc)).sqrt() * (nu / 2.0).tan().atan()
}

fn nu_to_F(nu: f64, ecc: f64) -> f64 {
    2.0 * ((ecc - 1.0) / (ecc + 1.0)).sqrt() * (nu / 2.0).tan().atanh()
}

fn E_to_nu(E: f64, ecc: f64) -> f64 {
    2.0 * ((1.0 + ecc) / (1.0 - ecc)).sqrt() * (E / 2.0).tan().atan()
}

fn F_to_nu(F: f64, ecc: f64) -> f64 {
    2.0 * ((ecc + 1.0) / (ecc - 1.0)).sqrt() * (F / 2.0).tanh().atan()
}

fn M_to_E(M: f64, ecc: f64) -> f64 {
    let E0 = if -PI < M && M < 0.0 || M > PI {
        M - ecc
    } else {
        M + ecc
    };
    newton_elliptic(E0, M, ecc, 1.48e-08, 50)
}

fn M_to_F(M: f64, ecc: f64) -> f64 {
    let F0 = (M / ecc).asinh();
    newton_hyperbolic(F0, M, ecc, 1.48e-08, 100)
}

fn M_to_D(M: f64) -> f64 {
    let B = 3.0 * M / 2.0;
    let A = (B + (1.0 + B * B).sqrt()).powf(2.0 / 3.0);
    2.0 * A * B / (1.0 + A + A * A)
}

fn E_to_M(E: f64, ecc: f64) -> f64 {
    E - ecc * E.sin()
}

fn F_to_M(F: f64, ecc: f64) -> f64 {
    ecc * F.sinh() - F
}

fn D_to_M(D: f64) -> f64 {
    D + D.powi(3) / 3.0
}

fn fp_angle(nu: f64, ecc: f64) -> f64 {
    (ecc * nu.sin()).atan2(1.0 + ecc * nu.cos())
}

fn kepler_equation_near_parabolic(D: f64, M: f64, ecc: f64) -> f64 {
    D_to_M_near_parabolic(D, ecc) - M
}

fn kepler_equation_prime_near_parabolic(D: f64, _M: f64, ecc: f64) -> f64 {
    let x = (ecc - 1.0) / (ecc + 1.0) * (D * D);
    assert!(x.abs() < 1.0);
    let S = dS_x_alt(ecc, x, 1e-12);
    (2.0 / (1.0 + ecc)).sqrt() + (2.0 / (1.0 + ecc).powi(3)).sqrt() * (D * D) * S
}

fn S_x(ecc: f64, x: f64, atol: f64) -> f64 {
    assert!(x.abs() < 1.0);
    let mut S = 0.0;
    let mut k = 0;
    loop {
        let S_old = S;
        S += (ecc - 1.0 / (2.0 * k as f64 + 3.0)) * x.powi(k);
        k += 1;
        if (S - S_old).abs() < atol {
            return S;
        }
    }
}

fn dS_x_alt(ecc: f64, x: f64, atol: f64) -> f64 {
    assert!(x.abs() < 1.0);
    let mut S = 0.0;
    let mut k = 0;
    loop {
        let S_old = S;
        S += (ecc - 1.0 / (2.0 * k as f64 + 3.0)) * (2.0 * k as f64 + 3.0) * x.powi(k);
        k += 1;
        if (S - S_old).abs() < atol {
            return S;
        }
    }
}

fn d2S_x_alt(ecc: f64, x: f64, atol: f64) -> f64 {
    assert!(x.abs() < 1.0);
    let mut S = 0.0;
    let mut k = 0;
    loop {
        let S_old = S;
        S += (ecc - 1.0 / (2.0 * k as f64 + 3.0))
            * (2.0 * k as f64 + 3.0)
            * (2.0 * k as f64 + 2.0)
            * x.powi(k);
        k += 1;
        if (S - S_old).abs() < atol {
            return S;
        }
    }
}

fn D_to_M_near_parabolic(D: f64, ecc: f64) -> f64 {
    let x = (ecc - 1.0) / (ecc + 1.0) * (D * D);
    assert!(x.abs() < 1.0);
    let S = S_x(ecc, x, 1e-12);
    (2.0 / (1.0 + ecc)).sqrt() * D + (2.0 / (1.0 + ecc).powi(3)).sqrt() * (D * D * D) * S
}

fn M_to_D_near_parabolic(M: f64, ecc: f64, tol: f64, maxiter: usize) -> f64 {
    let mut D0 = M_to_D(M);

    for _ in 0..maxiter {
        let fval = kepler_equation_near_parabolic(D0, M, ecc);
        let fder = kepler_equation_prime_near_parabolic(D0, M, ecc);

        let newton_step = fval / fder;
        let D = D0 - newton_step;
        if (D - D0).abs() < tol {
            return D;
        }

        D0 = D;
    }

    f64::NAN
}

fn delta_t_from_nu(nu: f64, ecc: f64, k: f64, q: f64, delta: f64) -> f64 {
    assert!(-PI <= nu && nu < PI);
    if ecc < 1.0 - delta {
        let E = nu_to_E(nu, ecc);
        let M = E_to_M(E, ecc);
        let n = (k * (1.0 - ecc).powi(3) / q.powi(3)).sqrt();
        M / n
    } else if 1.0 - delta <= ecc && ecc < 1.0 {
        let E = nu_to_E(nu, ecc);
        if delta <= 1.0 - ecc * E.cos() {
            let M = E_to_M(E, ecc);
            let n = (k * (1.0 - ecc).powi(3) / q.powi(3)).sqrt();
            M / n
        } else {
            let D = nu_to_D(nu);
            let M = D_to_M_near_parabolic(D, ecc);
            let n = (k / (2.0 * q.powi(3))).sqrt();
            M / n
        }
    } else if ecc == 1.0 {
        let D = nu_to_D(nu);
        let M = D_to_M(D);
        let n = (k / (2.0 * q.powi(3))).sqrt();
        M / n
    } else if 1.0 + ecc * nu.cos() < 0.0 {
        f64::NAN
    } else if 1.0 < ecc && ecc <= 1.0 + delta {
        let F = nu_to_F(nu, ecc);
        if delta <= ecc * F.cosh() - 1.0 {
            let M = F_to_M(F, ecc);
            let n = (k * (ecc - 1.0).powi(3) / q.powi(3)).sqrt();
            M / n
        } else {
            let D = nu_to_D(nu);
            let M = D_to_M_near_parabolic(D, ecc);
            let n = (k / (2.0 * q.powi(3))).sqrt();
            M / n
        }
    } else if 1.0 + delta < ecc {
        let F = nu_to_F(nu, ecc);
        let M = F_to_M(F, ecc);
        let n = (k * (ecc - 1.0).powi(3) / q.powi(3)).sqrt();
        M / n
    } else {
        panic!("Unreachable code");
    }
}

fn nu_from_delta_t(delta_t: f64, ecc: f64, k: f64, q: f64, delta: f64) -> f64 {
    if ecc < 1.0 - delta {
        let n = (k * (1.0 - ecc).powi(3) / q.powi(3)).sqrt();
        let M = n * delta_t;
        let E = M_to_E((M + PI) % (2.0 * PI) - PI, ecc);
        E_to_nu(E, ecc)
    } else if 1.0 - delta <= ecc && ecc < 1.0 {
        let E_delta = ((1.0 - delta) / ecc).acos();
        let n = (k * (1.0 - ecc).powi(3) / q.powi(3)).sqrt();
        let M = n * delta_t;
        if E_to_M(E_delta, ecc).abs() <= M.abs() {
            let E = M_to_E((M + PI) % (2.0 * PI) - PI, ecc);
            E_to_nu(E, ecc)
        } else {
            let n = (k / (2.0 * q.powi(3))).sqrt();
            let M = n * delta_t;
            let D = M_to_D_near_parabolic(M, ecc, 1.48e-8, 50);
            D_to_nu(D)
        }
    } else if ecc == 1.0 {
        let n = (k / (2.0 * q.powi(3))).sqrt();
        let M = n * delta_t;
        let D = M_to_D(M);
        D_to_nu(D)
    } else if 1.0 < ecc && ecc <= 1.0 + delta {
        let F_delta = ((1.0 + delta) / ecc).acosh();
        let n = (k * (ecc - 1.0).powi(3) / q.powi(3)).sqrt();
        let M = n * delta_t;
        if F_to_M(F_delta, ecc).abs() <= M.abs() {
            let F = M_to_F(M, ecc);
            F_to_nu(F, ecc)
        } else {
            let n = (k / (2.0 * q.powi(3))).sqrt();
            let M = n * delta_t;
            let D = M_to_D_near_parabolic(M, ecc, 1.48e-8, 50);
            D_to_nu(D)
        }
    } else {
        let n = (k * (ecc - 1.0).powi(3) / q.powi(3)).sqrt();
        let M = n * delta_t;
        let F = M_to_F(M, ecc);
        F_to_nu(F, ecc)
    }
}

fn farnocchia_coe(
    k: f64,
    p: f64,
    ecc: f64,
    inc: f64,
    raan: f64,
    argp: f64,
    nu: f64,
    tof: f64,
) -> f64 {
    let q = p / (1.0 + ecc);

    let delta_t0 = delta_t_from_nu(nu, ecc, k, q, 1e-2);
    let delta_t = delta_t0 + tof;

    nu_from_delta_t(delta_t, ecc, k, q, 1e-2)
}

struct Orbit {
    k: f64,
    p: f64,
    ecc: f64,
    inc_deg: f64,
    raan_deg: f64,
    argp_deg: f64,
    nu_deg: f64,
}

#[cfg(test)]
mod test {
    use super::*;

    use rstest::*;

    #[fixture]
    fn orbit() -> Orbit {
        Orbit {
            k: 398600.4418,
            p: 6780.8472106,
            ecc: 0.00130547,
            inc_deg: 51.6012092,
            raan_deg: 198.37949974,
            argp_deg: 39.26289661,
            nu_deg: 46.59580468,
        }
    }

    #[rstest]
    fn test_noop_farnocchia_coe(orbit: Orbit) {
        let tof = 0.0;
        let expected_nu: f64 = orbit.nu_deg.to_radians();

        let nu = farnocchia_coe(
            orbit.k,
            orbit.p,
            orbit.ecc,
            orbit.inc_deg.to_radians(),
            orbit.raan_deg.to_radians(),
            orbit.argp_deg.to_radians(),
            orbit.nu_deg.to_radians(),
            tof,
        );

        assert!((nu - expected_nu).abs() < 1e-15);
    }

    #[rstest]
    fn test_short_farnocchia_coe(orbit: Orbit) {
        let tof = 100.0;
        let expected_nu: f64 = 0.9265094210290502; // Computed with poliastro
        //       computed_nu = 0.9265619416308579  // TODO: Can it be closer?

        let nu = farnocchia_coe(
            orbit.k,
            orbit.p,
            orbit.ecc,
            orbit.inc_deg.to_radians(),
            orbit.raan_deg.to_radians(),
            orbit.argp_deg.to_radians(),
            orbit.nu_deg.to_radians(),
            tof,
        );

        println!("Expected nu: {}", expected_nu);
        println!("Computed nu: {}", nu);
        assert!((nu - expected_nu).abs() < 1e-4);
    }

    #[rstest]
    fn test_period_farnocchia_coe(orbit: Orbit) {
        let tof = 5556.969701163016; // Computed with poliastro
        let expected_nu: f64 = orbit.nu_deg.to_radians();
        //       computed_nu = 0.8132502093265356  // TODO: Can it be closer?

        let nu = farnocchia_coe(
            orbit.k,
            orbit.p,
            orbit.ecc,
            orbit.inc_deg.to_radians(),
            orbit.raan_deg.to_radians(),
            orbit.argp_deg.to_radians(),
            orbit.nu_deg.to_radians(),
            tof,
        );

        println!("Expected nu: {}", expected_nu);
        println!("Computed nu: {}", nu);
        assert!((nu - expected_nu).abs() < 1e-10);
    }

    #[rstest]
    fn test_long_farnocchia_coe(orbit: Orbit) {
        let tof = 20000.;
        let expected_nu: f64 = -1.7102617293760711; // Computed with poliastro
        //       computed_nu = -1.7113142537561092  // TODO: Can it be closer?

        let nu = farnocchia_coe(
            orbit.k,
            orbit.p,
            orbit.ecc,
            orbit.inc_deg.to_radians(),
            orbit.raan_deg.to_radians(),
            orbit.argp_deg.to_radians(),
            orbit.nu_deg.to_radians(),
            tof,
        );

        println!("Expected nu: {}", expected_nu);
        println!("Computed nu: {}", nu);
        assert!((nu - expected_nu).abs() < 1e-2);
    }
}
