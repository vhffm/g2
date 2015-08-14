"""
Impact Simulation Helpers.
"""

import numpy as np
import constants as C
import scipy as sp
import io_helpers as ioh
import pandas as pd


def calibrate_craters_n83(dfc_all, blk_all, tlo_all, thi_all, pid_target=2):
    """
    Calibrate Collision Rate to Neukum (1983) Cratering Fits.
    Time Range: tlo_all =< time =< thi_all
    
    Cf. Neukum+ 2001, Fig. 10 / Eqn. 5
    http://link.springer.com/article/10.1023/A:1011989004263

    @param: dfc_all - List of Collision Logs [List o/ Pandas Dataframe]
    @param: tlo_all - List of Lower Time Limit to Calibrate Against (Years) 
                      [Numpy Float Array]
    @param: thi_all - List of Upper Time Limit to Calibrate Against (Years)
                      [Numpy Float Array]
    @param: blk_all - List of Blacklisted Particles [List o/ Numpy Arrays]
    @param: pid_target - Particle ID of Earth (2 or 3) [Integer]
    @return: crc_all - Correction Factors [Numpy Float Array]
    """

    crc_all = []
    for idfc, dfc in enumerate(dfc_all):
        dfc = dfc[dfc.pidi==pid_target]

        # Only Consider Subset in Time
        if not np.isnan(tlo_all[idfc]):
            # print "Lower Limit - %i - %.2e" % (idfc, tlo_all[idfc])
            dfc = dfc[dfc.time >= tlo_all[idfc]]

        if not np.isnan(thi_all[idfc]):
            # print "Upper Limit - %i - %.2e" % (idfc, thi_all[idfc])
            dfc = dfc[dfc.time <= thi_all[idfc]]

        # Time
        time = dfc[~dfc.pidj.isin(blk_all[idfc]) & \
                   (dfc.pidi==pid_target)].time[::-1]
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


def neukum_production_function(Dmin=0.01, Dmax=300.0):
    """
    Neukum Crater Production Distribution. Based on Neukum 1983.
    Cf. Neukum+ 2001, Eqn. 2, Eqn. 3, Table 1
    http://link.springer.com/article/10.1023/A:1011989004263

    *** WARNING ***
    Valid in the range 0.01 - 300.0 km diameter. Extrapolate carefully.
    Cummulative crater count in calibrate_craters_n83() valid D > 1km only!
    For details, see Neukum+ 2001, Page 15, Text Below Eqn. 5

    @param: Dmin - Lower End of Crater Mass Function [Numpy Float]
    @param: Dmax - Upper End of Crater Mass Function [Numpy Float]
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
    D = np.logspace(np.log10(Dmin), np.log10(Dmax), 4096) # km

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
    return M, L


def scale_production_function(crc_all, Dmin=1.0, Dmax=300.0):
    """
    Computes Number Scaling from Production Functions.

    The cummulative crater count we scale to is valid only for D > 1 km.
    Cf. calibrate_craters_n83() and neukum_production_function() for comments.

    @param: crc_all - N of Real Impacts per Sim Impact [Numpy Float Array]
    @param: Dmin - Lower End of Crater Mass Function [Numpy Float]
    @param: Dmax - Upper End of Crater Mass Function [Numpy Float]
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
    _, dNdD, _, D = neukum_production_function(Dmin, Dmax)
    actual = sp.integrate.simps(dNdD, x=D)

    # We calibrate our counted impacts to the expected impacts.
    nscale = counted / actual

    # Using a typical impactor mass for a given diameter D, we compute the 
    # the total mass delivered over these impacts. This sets the mass of our 
    # simulation particles. Impact velocity and angle are median values
    # taken from Test Particle simulations.
    M, _ = impactor_mass(D, 12.0, 18.5 * C.d2r)
    mscale = sp.integrate.simps((dNdD*M), x=D)
    mscale *= nscale
    mscale *= C.Smoon * (C.Aearth/C.Amoon)

    # Return Scaling
    return nscale, mscale


def load(run):
    """
    Load a given simulation. Beware of hardcoded paths.

    @param: run - Name of run to load [String]
    @return: dfc - Collisions [Pandas Dataframe]
    @return: dfe - Ejections  [Pandas Dataframe]
    @return: dfo - Coordinate Outputs [Pandas Dataframe]
    @return: blacklist - Blacklisted Partilces [Numpy Integer Array]
    """

    # Datadir
    datadir = '/zbox/data/volker/Impacts'

    # Which Simulation?
    if run == 'morby':
        basedir = "%s/Morby" % datadir
        blckdir = basedir
        run_name = 'impacts_morby'
    elif run == 'morby_hd':
        basedir = "%s/Morby_HighDispersion" % datadir
        blckdir = basedir
        run_name = 'impacts_morby_hd'
    elif run == 'solar2':
        basedir = "%s/Solar2" % datadir
        blckdir = basedir
        run_name = 'impacts_solar2_long'
    elif run == 'solar2_hd':
        basedir = "%s/Solar2_HighDispersion" % datadir
        blckdir = basedir
        run_name = 'impacts_solar2_hd'
    elif run == 'blowup':
        basedir = "%s/Blowup" % datadir
        blckdir = "%s/Morby" % datadir
        run_name = 'impacts_morby'
    elif run == 'blowup2':
        basedir = "%s/Blowup2" % datadir
        blckdir = "%s/Morby" % datadir
        run_name = 'impacts_morby'
    elif run == 'blowup2_hd':
        basedir = "%s/Blowup2_HighDispersion" % datadir
        blckdir = "%s/Morby_HighDispersion" % datadir
        run_name = 'impacts_morby_hd'
    else:
        raise Exception("Invalid Run %s" % run)

    # Coordinate Output Steps & Number of Runs
    if run in [ 'morby', 'morby_hd' ]:
        nsteps = np.array([0,1,2,3,4,5], dtype=np.int64) * 1e9
        nrun_hi = 16
    elif run in [ 'solar2' ]:
        nsteps = np.array([0,1,2,3,4,5,6,7], dtype=np.int64) * 1e9
        nrun_hi = 16
    elif run in [ 'solar2_hd' ]:
        nsteps = np.array([0,1,2,3,4,5,6,7,8,9,10], dtype=np.int64) * 1e9
        nrun_hi = 16
    elif run in [ 'blowup', 'blowup2', 'blowup2_hd' ]:
        nsteps = np.array([5,6,7,8,9,10], dtype=np.int64) * 1e9
        nrun_hi = 16
    else:
        raise Exception("Invalid Run %s" % run)
        
    # Debug
    print "// Loading %s" % run

    # Number of Runs?
    # 1..16 are default
    #    17 is 0.5 - 2.0 AU
    nruns = range(1,nrun_hi+1)

    fnames_c = []
    fnames_e = []
    fnames_o = []
    for nrun in nruns:
        fnames_c.append("%s/%02d/Collisions_%s_%02d.dat" % \
            (basedir, nrun, run_name, nrun))
        fnames_e.append("%s/%02d/Ejections_%s_%02d.dat" % \
            (basedir, nrun, run_name, nrun))
        for _, nstep in enumerate(nsteps):
            fnames_o.append("%s/%02d/Out_%s_%02d_%012d.dat" % \
                (basedir, nrun, run_name, nrun, nstep))

    dfc = ioh.read_collisions_and_stack(fnames_c, return_xyz=True)
    dfe = ioh.read_ejections_and_stack(fnames_e)
    dfo = ioh.read_output_and_stack(fnames_o, frame='heliocentric')

    dfc.sort(columns=["time"], inplace=True)
    dfe.sort(columns=["time"], inplace=True)
    
    # Chop, Rewind?
    if run in [ 'blowup', 'blowup2', 'blowup2_hd' ]:
        tblowup = 5.0e9 * 36.0/365.25
        dfc = dfc[dfc.time>tblowup]
        dfe = dfe[dfe.time>tblowup]
        # dfc["time"] -= tblowup
        # dfe["time"] -= tblowup

    # Units
    dfc["theta"] *= C.r2d

    # Blacklist
    blacklist = np.load("%s/Blacklist_1Rhill.npz" % blckdir)["blacklist"]
    
    # Time Filter
    # dfc = dfc[dfc.time>1.0e7]
    # dfe = dfe[dfe.time>1.0e7]

    # Return
    return dfc, dfe, dfo, blacklist


def load_all():
    """
    Load all simulations. Also join Morby/Blowup2, Morby_HD/Blowup2_HD.
    Hardcoding is bad, but so are 10 monstrous cells in IPython just to load.

    @return: dfc_all - Dataframes w/ Collisions [List of Dataframes]
    @return: dfe_all - Dataframes w/ Ejections  [List of Dataframes]
    @return: dfo_all - Dataframes w/ Coordinate Outputs [List of Dataframes]
    @return: blk_all - Blacklisted Particle IDs [List of Numpy Arrays]
    @return: tag_all - Simulation Tags [List of Strings]
    @return: tlo_all - Lower Time Bound for Calibration [Numpy Float Array]
    @return: thi_all - Upper Time Bound for Calibration [Numpy Float Array]
    """
    
    # Load
    dfc_solar2, dfe_solar2, dfo_solar2, blacklist_solar2 = \
        load('solar2')
    dfc_solar2_hd, dfe_solar2_hd, dfo_solar2_hd, blacklist_solar2_hd = \
        load('solar2_hd')
    dfc_morby, dfe_morby, dfo_morby, blacklist_morby = \
        load('morby')
    dfc_morby_hd, dfe_morby_hd, dfo_morby_hd, blacklist_morby_hd = \
        load('morby_hd')
    dfc_blowup, dfe_blowup, dfo_blowup, _ = \
        load('blowup')
    dfc_blowup2, dfe_blowup2, dfo_blowup2, _ = \
        load('blowup2')
    dfc_blowup2_hd, dfe_blowup2_hd, dfo_blowup2_hd, _ = \
        load('blowup2_hd')
    
    # Blowup Time
    tblowup = 5.0e9 * 36.0/365.25
    
    # Calibration Target (Earth)
    pid_target = 2
    
    # Master List
    dfo_all = [ dfo_solar2, dfo_solar2_hd, \
                dfo_morby, dfo_morby_hd, \
                dfo_blowup, dfo_blowup2, dfo_blowup2_hd ]
    dfc_all = [ dfc_solar2, dfc_solar2_hd, \
                dfc_morby, dfc_morby_hd, \
                dfc_blowup, dfc_blowup2, dfc_blowup2_hd ]
    dfe_all = [ dfe_solar2, dfe_solar2_hd, \
                dfe_morby, dfe_morby_hd, \
                dfe_blowup, dfe_blowup2, dfe_blowup2_hd ]
    blk_all = [ blacklist_solar2, blacklist_solar2_hd, \
                blacklist_morby, blacklist_morby_hd, \
                np.array([]), np.array([]), np.array([]) ]
    tag_all = [ "Solar2", "Solar2/HD", \
                "Morby", "Morby/HD", \
                "Blowup", "Blowup2", "Blowup2/HD" ]

    #
    # Time Range for Calibration
    #
    # NB:
    # For Solar2/Solar2_HD runs, the lower time limit for calibration is the 
    # time after the blow up of the model. ALthough there is no blow-up in
    # this model, we only have actual evidence for the time after this. So 
    # there's no point in calibrating to the time before.
    #
    # For Morby/Morby_HD runs, the lower time limit is a dummy values. All
    # calibrations for these are obtained from Blowup2/Blowup2_HD runs
    #
    # For the upper time limit, we chop off the last 100 Myr because our
    # cummulative impact count drops steeply. Courtesy of being logarithmic.
    #
    tlo_all = np.array([tblowup, tblowup, \
                        0.0, 0.0, \
                        tblowup, tblowup, tblowup ])
    thi_all = np.ones_like(tlo_all) * 9.0e9 * 36.0 / 365.25
    
    assert len(dfc_all) == len(dfe_all) == len(blk_all) == \
        len(tag_all) == len(tlo_all) == len(thi_all) == len(dfo_all)
    
    ###########################################################################

    # A - Tack Together Morby + Blowup2
    dfc1 = dfc_all[2].copy()
    dfc2 = dfc_all[5].copy()
    dfc1x = dfc1
    dfc2x = dfc2
    # dfc1x = dfc1[~dfc1.pidj.isin(dfc1.pidi==pid_target)] # morby
    # dfc2x = dfc2[~dfc2.pidj.isin(dfc2.pidi==pid_target)] # blowup2
    # dfc2x.time += tblowup
    dfcxx = pd.concat([dfc1x, dfc2x])
    dfoxx = pd.concat([dfo_all[2].copy(), dfo_all[5].copy()])
    dfoxx.drop_duplicates(subset=["pid", "time", "nstep"], inplace=True)
    dfoxx.reset_index(drop=True, inplace=True)

    # Append to Lists
    tlo_all = np.append(tlo_all, tblowup)
    thi_all = np.append(thi_all, np.nan)
    dfo_all.append(dfoxx)
    dfc_all.append(dfcxx)
    dfe_all.append(pd.DataFrame())
    blk_all.append(blacklist_morby)
    tag_all.append("Morby + Blowup2")

    ###########################################################################

    # B - Tack Together Morby_HD + Blowup2_HD
    dfc1 = dfc_all[3].copy()
    dfc2 = dfc_all[6].copy()
    dfc1x = dfc1
    dfc2x = dfc2
    # dfc1x = dfc1[~dfc1.pidj.isin(dfc1.pidi==pid_target)] # morby_hd
    # dfc2x = dfc2[~dfc2.pidj.isin(dfc2.pidi==pid_target)] # blowup2_hd
    # dfc2x.time += tblowup
    dfcxx = pd.concat([dfc1x, dfc2x])
    dfoxx = pd.concat([dfo_all[3].copy(), dfo_all[6].copy()])
    dfoxx.drop_duplicates(subset=["pid", "time", "nstep"], inplace=True)
    dfoxx.reset_index(drop=True, inplace=True)

    # Append to Lists
    tlo_all = np.append(tlo_all, tblowup)
    thi_all = np.append(thi_all, np.nan)
    dfo_all.append(dfoxx)
    dfc_all.append(dfcxx)
    dfe_all.append(pd.DataFrame())
    blk_all.append(blacklist_morby)
    tag_all.append("Morby/HD + Blowup2/HD")

    ###########################################################################

    return dfc_all, dfe_all, dfo_all, blk_all, tag_all, tlo_all, thi_all


def calibrate_all(dfc_all, blk_all, tlo_all, thi_all, Dmin=1.0, Dmax=300.0):
    """
    Calibrate simulations against Neukum+ 1983 cummulative crater count.
    Note that Morby(_HD) uses Blowup2(_HD) calibration.

    The number count correction (crc_all) gives the number of actual particles
    (in some diameter range) that one simulation particle corresponds to.

    The mass correction (mscale) gives the total mass each simulation particle
    corresponds to.

    Both assume a power-law-ish size-frequency distribution of particles. 
    Cf. calibrate_craters_n83() and scale_production_function() for refs.

    The cummulative crater count is valid only for D > 1 km.
    Cf. calibrate_craters_n83() and neukum_production_function() for comments.

    @param: dfc_all - Dataframes w/ Collisions [List of Dataframes]
    @param: blk_all - Blacklisted Particle IDs [List of Numpy Arrays]
    @param: tlo_all - Lower Time Bound for Calibration [Numpy Float Array]
    @param: thi_all - Upper Time Bound for Calibration [Numpy Float Array]
    @param: Dmin - Lower End of Crater Mass Function [Numpy Float]
    @param: Dmax - Upper End of Crater Mass Function [Numpy Float]
    @return: crc_all - Number Count Correction Factors [Numpy Float Array]
    @return: mscale  - Mass Correction Factors [Numpy Float Array]
    """
    
    # Compute Offsets
    crc_all = calibrate_craters_n83(dfc_all, blk_all, tlo_all, thi_all)

    # Blowup2 => Morby
    crc_all[2] = crc_all[5]

    # Blowup2_HD => Morby_HD
    crc_all[3] = crc_all[6]

    # Blowup2 => XXX
    # crc_all = np.append(crc_all, crc_all[2])
    # crc_all = np.append(crc_all, crc_all[3])

    ###########################################################################

    # Compute Mass Scale
    _, mscale = scale_production_function(crc_all, Dmin, Dmax)
    
    ###########################################################################
    
    return crc_all, mscale
