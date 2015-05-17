"""
Impact Simulation Helpers.
"""

import numpy as np
import constants as C
import scipy as sp


def calibrate_craters_n83(dfc_all, t00_all, blk_all, pid_target=2):
    """
    Calibrate Collision Rate to Neukum (1983) Cratering Fits.
    Cf. Neukum+ 2001, Fig. 10 / Eqn. 5
    http://link.springer.com/article/10.1023/A:1011989004263

    @param: dfc_all - List of Collision Logs [List o/ Pandas Dataframe]
    @param: t00_all - List of Initial Times (Years) [Numpy Float Array]
    @param: blk_all - List of Blacklisted Particles [List o/ Numpy Arrays]
    @param: pid_target - Particle ID of Earth (2 or 3) [Integer]
    @return: crc_all - Correction Factors [Numpy Float Array]
    """

    crc_all = []
    for idfc, dfc in enumerate(dfc_all):
        dfc = dfc[dfc.pidi==pid_target]

        # Time
        time = dfc[~dfc.pidj.isin(blk_all[idfc]) & \
                   (dfc.pidi==pid_target)].time[::-1]
        time += t00_all[idfc] # Beginning of Sim
        time /= 1.0e9         # Yr => Gyr
        time += 0.115         # Gyr
        time += C.t0ss        # Beginning of Solar System
        time = np.asarray(time)

        # Cummulative Collision Ccount
        coll = np.cumsum(np.ones(len(dfc[~dfc.pidj.isin(blk_all[idfc]) & \
                                         (dfc.pidi==pid_target)])))
        coll *= C.Amoon / C.Aearth
        coll /= C.Smoon

        # Neukum Count
        time_n83 = np.linspace(time[0], time[-1], 512)
        coll_n83 = 5.44e-14 * ( np.exp(6.93 * (-time_n83)) - 1 ) + \
            8.3e-4 * (-time_n83)
        
        # Interpolate Collision to Same Grid
        coll_interp = np.interp(time_n83[::-1], time[::-1], coll[::-1])
        coll_interp = coll_interp[::-1]
        
        # Output Correction
        crc_all.append(np.median(coll_n83/coll_interp))

    # Array
    crc_all = np.asarray(crc_all)

    # Return
    return crc_all


def neukum_production_function():
    """
    Neukum Crater Production Distribution. Based on Neukum 1983.
    Cf. Neukum+ 2001, Eqn. 2, Eqn. 3, Table 1
    http://link.springer.com/article/10.1023/A:1011989004263

    @return: N - Numnber of Craters @ Diameter [Numpy Float Array]
    @return: dNdR - Diff. Number of Craters @ Diameter [Numpy Float Array]
    @return: R - R-Distribution [Numpy Float Array]
    @return: D - Diameters [Numpy Array]
    """

    a0 = -3.0876; a1 = -3.557528; a2 = 0.781027; a3 = 1.021512
    a4 = -0.156012; a5 = -0.444058; a6 = 0.019977; a7 = 0.086850
    a8 = -0.005874; a9 = -0.006809; a10 = 8.25e-4; a11 = 5.54e-5

    a = np.array([a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11])
    # D = np.logspace(-2, 3, 4096) # km
    D = np.logspace(-2, np.log10(300.0), 4096) # km

    # N(D) -- Absolute Distribution
    logN = np.zeros_like(D)
    for ii in range(len(a)):
        logN += a[ii] * np.log10(D)**ii
    N = 10.0**logN

    # dN/dD -- Differential Distribution
    logNx = np.zeros_like(D)
    for ii in range(1, len(a)):
        logNx += a[ii] * np.log10(D)**(ii-1)
    dNdD = (N/D) * logNx
    dNdD = np.abs(dNdD)

    # R-Distribution
    R = dNdD * D**3.0

    # Return
    return N, dNdD, R, D


def impactor_mass(D, v_i, theta, rho_t=3.34, rho_i=2.00, grav=1.622):
    """
    Compute Size & Mass of Impactor for Crater Diameter D.
    Cf. Collins+ 2005, Eqn. 21
    http://adsabs.harvard.edu/abs/2005M&PS...40..817C
    http://impact.ese.ic.ac.uk/ImpactEffects/effects.pdf

    Default target is Moon. Impactor density is mean asteroid density.
    Test particle simulations indicate typical values for v_i ~ 12 km/s and
    theta ~ 18.5 Degree.

    @param: D - Crater diameter (km) [Float or Numpy Float Array]
    @param: vi - Impact velocity (km/s) [Float]
    @param: theta - Impact angle (radian) [Float]
    @param: rho_t - Target density (g/cc) [Float]
    @param: rho_i - Impactor density (g/cc) [Float]
    @param: grav - Target surface gravity (m/s2) [Float]
    @return: M - Impactor mass (kg) [Float or Numpy Float Array]
    """

    # Fix Units
    grav /= 1000.0 # m/s2 => km/s2

    # Compute Impactor Diameter
    L = D
    L /= 1.161 * (rho_i/rho_t)**(1.0/3.0) * \
         v_i**0.44 * \
         grav**(-0.22) * \
         np.sin(theta)**(1.0/3.0)
    L = L**(1.0/0.78)

    # Compute (Spherical) Impactor Mass
    M = 4.0/3.0 * np.pi * (L/2.0)**3.0 * rho_i/1000.0*(1000.0*100.0)**(3.0)

    # Return [Float or Numpy Float Array; Same as D]
    return M


def scale_production_function(crc_all):
    """
    Computes Number Scaling from Production Functions.

    @param: crc_all - N of Real Impacts per Sim Impact [Numpy Float Array]
    @return: nscale - Resulting Production Function Scaling [Numpy Float Array]
    @return: mscale - Resulting Mass Function Scaling [Numpy Float Array]

    The effective mass of a given Simulation Particle is "mscale".
    """

    # Each Earth impact in the simulation corresponds to crc_all actual
    # impacts. We scale this to the Moon's cross section (Amoon/Aearth), 
    # and divide by Moon's surface area to get the impacts / km^2.
    counted = crc_all * (C.Amoon/C.Aearth) / C.Smoon

    # We integrate the differential Production Function to obtain the number
    # of expected impacts in the relevant mass ranges (1km < D < 300km).
    _, dNdD, _, D = neukum_production_function()
    actual = sp.integrate.simps(dNdD[D>1.0], x=D[D>1.0])

    # We calibrate our counted impacts to the expected impacts.
    nscale = counted / actual

    # Using a typical impactor mass for a given diameter D, we compute the 
    # the total mass delivered over these impacts. This sets the mass of our 
    # simulation particles. Impact velocity and angle are median values
    # taken from Test Particle simulations.
    M = impactor_mass(D, 12.0, 18.5 * C.d2r)
    mscale = sp.integrate.simps((dNdD*M)[D>1.0], x=D[D>1.0])
    mscale *= nscale
    mscale *= C.Smoon * (C.Aearth/C.Amoon)

    # Return Scaling
    return nscale, mscale
