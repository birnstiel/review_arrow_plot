#!/usr/bin/env python

import numpy as np
import scipy.constants as sc
pi = np.pi                   # PI
k_b = sc.k * 1e7                # Boltzmann constant in erg/K
m_p = sc.proton_mass * 1e3      # proton mass in g
Grav = sc.G * 1e3                # gravitational constant in cm^3 g^-1 s^-2
AU = sc.au * 1e2               # astronomical unit in cm
year = sc.Julian_year          # year in s
mu = 2.3e0                   # mean molecular mass in proton masses
M_sun = 1.9891e+33              # mass of the sun in g
sigma_H2 = 2e-15                    # H2 collisional cross section


def dlydlx(x, R):
    """
    calculates the log-derivative

     dlog(y)
    -------- = dlydlx(x,y)
     dlog(x)
    """
    from numpy import zeros, shape, interp, log10
    #
    # define the interpolation function (one for each row)
    #
    r = zeros(shape(R))
    if len(shape(R)) > 1:
        for i, row in enumerate(R):
            def R_int(x_int):
                return 10**(interp(log10(x_int), log10(x), log10(row)))
            h = x / 100.
            r[i] = x / row * (R_int(x + h) - R_int(x - h)) / (2. * h)
    else:
        def R_int(x_int):
            return 10**(interp(log10(x_int), log10(x), log10(R)))
        h = x / 100.
        r = x / R*(R_int(x + h) - R_int(x - h)) / (2. * h)
    return r
