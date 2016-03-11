"""
Helpers for Isotope Tracking.
"""

import numpy as np
import pandas as pd
import other_helpers as oh
import formation_helpers as fh


def assign_ti_concentration(a, model=1):
    """
    Assign Titanium Concentration.
    Concentration is Parts per Million (PPM).

    @param:  a       - Semi-Major Axis        (AU)  [Float]
    @param:  model   - Model Numner (1,2,3)         [Integer]
    @return: ti_conc - Titanium Concentration (PPM) [Float]
    """

    if not model in [ 1, 2, 3 ]:
        raise Exception('Invalid Model.')

    if model == 1:
        # @todo: Check w/ Maria. This is >2.8 AU in her mail, 
        #        but then 2.7 to 2.8 range is undefined
        if a >= 2.7:
            ti_conc = 450.0
        elif a >= 2.5 and a < 2.7:
            ti_conc = 900.0
        elif a >= 2.0 and a < 2.5:
            ti_conc = 600.0
        elif a < 2.0:
            ti_conc = 810.0

    elif model == 2:
        if a < 2.7:
            m, n = oh.mkline(0.5, 810.0, 2.7, 900.0)
            ti_conc = m * a + n
        elif a >= 2.7 and a < 2.8:
            m, n = oh.mkline(2.7, 900.0, 2.8, 450.0)
            ti_conc = m * a + n
        elif a >= 2.8:
            ti_conc = 450.0

    elif model == 3:
        m, n = oh.mkline(0.5, 810.0, 4.0, 270.0)
        ti_conc = m * a + n

    # Rescale to Fraction
    # ti_conc /= 1.0e6

    # Return
    return ti_conc


def assign_e50ti(a, model=1):
    """
    Assign Titanium Isotope Ratio e(50Ti/47Ti).

    @param:  a      - Semi-Major Axis (AU)    [Float]
    @param:  model  - Model Numner    (1,2,3) [Integer]
    @return: e50ti  - Isotope Ratio   (-)     [Float]
    """

    if not model in [ 1, 2, 3 ]:
        raise Exception('Invalid Model.')

    if model == 1:
        # @todo: Check w/ Maria. This is >2.8 AU in her mail, 
        #        but then 2.7 to 2.8 range is undefined
        if a >= 2.7:
            e50ti = 1.7
        elif a >= 2.5 and a < 2.7:
            e50ti = 3.3
        elif a >= 2.0 and a < 2.5:
            e50ti = -0.52
        elif a < 2.0:
            e50ti = -1.3

    elif model == 2:
        if a < 2.7:
            m, n = oh.mkline(0.5, -2.0, 2.7, 3.3)
            e50ti = m * a + n
        elif a >= 2.7 and a < 2.8:
            m, n = oh.mkline(2.7, 3.3, 2.8, 1.79)
            e50ti = m * a + n
        elif a >= 2.8:
            e50ti = 1.79

    elif model == 3:
        m, n = oh.mkline(0.5, -2.0, 4.0, 10.0)
        e50ti = m * a + n

    # Return
    return e50ti


def compute_isotope_fractions(dfo, dfo_t0, dfc, showstep=False):
    """
    Compute Isotope Fractions.
    Build Source List for Output, Compute Fractions from Precursors.

    Basically a clone of the compute_wmf() defined above.

    @param:  dfo    - Coordinate Output @ Time               [Pandas Dataframe]
    @param:  dfo_t0 - Outout @ Initial Time                  [Pandas Dataframe]
    @param:  dfc    - Collision List                         [Pandas Dataframe]
    @return: dfo    - Coordinate Output @ Time w/ WMF Fields [Pandas Dataframe]
    """

    ti_conc_01_ppm = np.ones(len(dfo)) * np.nan
    ti_conc_02_ppm = np.ones(len(dfo)) * np.nan
    ti_conc_03_ppm = np.ones(len(dfo)) * np.nan
    
    e50ti_01 = np.ones(len(dfo)) * np.nan
    e50ti_02 = np.ones(len(dfo)) * np.nan
    e50ti_03 = np.ones(len(dfo)) * np.nan

    dfc_now = dfc[dfc.time<=dfo.time.iloc[0]]
    
    # The Murder Loop
    # @todo - Accelerate? Rewrite Source List Construction in Fortran?
    for ii, dfo_row in enumerate(dfo.iterrows()):
        # Info
        if showstep:
            if ii % int(len(dfo)/8) == 0:
                print "%i/%i" % (ii,len(dfo))
        
        # Extract Series from One-Row Dataframe
        dfo_loc = dfo_row[1]
        
        # Identify Source Particles
        sources = \
            fh.return_sources(int(dfo_loc.pid), \
                           dfc_now[dfc_now.ifname == int(dfo_loc.ifname)])
        dfo_sources = \
            dfo_t0[(dfo_t0.pid.isin(sources)) & \
                   (dfo_t0.ifname == int(dfo_loc.ifname))].copy()
        
        # Imprint Titanium Concentration
        dfo_sources.loc[:,'ti_conc_01_ppm'] = \
            pd.Series(dfo_sources.a.apply(assign_ti_concentration, model=1), \
                      index=dfo_sources.index)
        dfo_sources.loc[:,'ti_conc_02_ppm'] = \
            pd.Series(dfo_sources.a.apply(assign_ti_concentration, model=2), \
                      index=dfo_sources.index)
        dfo_sources.loc[:,'ti_conc_03_ppm'] = \
            pd.Series(dfo_sources.a.apply(assign_ti_concentration, model=3), \
                      index=dfo_sources.index)
            
        # Imprint e(50Ti/47Ti)
        dfo_sources.loc[:,'e50ti_01'] = \
            pd.Series(dfo_sources.a.apply(assign_e50ti, model=1), \
                      index=dfo_sources.index)
        dfo_sources.loc[:,'e50ti_02'] = \
            pd.Series(dfo_sources.a.apply(assign_e50ti, model=2), \
                      index=dfo_sources.index)
        dfo_sources.loc[:,'e50ti_03'] = \
            pd.Series(dfo_sources.a.apply(assign_e50ti, model=3), \
                      index=dfo_sources.index)
            
        # Weighted Sum for Titanium Concentration
        ti_conc_01_ppm[ii] = np.sum(dfo_sources.mass * \
                                    dfo_sources.ti_conc_01_ppm) / \
                             np.sum(dfo_sources.mass)
        ti_conc_02_ppm[ii] = np.sum(dfo_sources.mass * \
                                    dfo_sources.ti_conc_02_ppm) / \
                             np.sum(dfo_sources.mass)
        ti_conc_03_ppm[ii] = np.sum(dfo_sources.mass * \
                                    dfo_sources.ti_conc_03_ppm) / \
                             np.sum(dfo_sources.mass)
            
        # Weighted Sum for e(50Ti/47Ti)
        # 0.07437 is the abundance of 47Ti
        # Cf. https://en.wikipedia.org/wiki/Isotopes_of_titanium
        f_01 = \
            (0.07437 * dfo_sources.ti_conc_01_ppm * dfo_sources.mass) / \
            np.sum(0.07437 * dfo_sources.ti_conc_01_ppm * dfo_sources.mass)
        f_02 = \
            (0.07437 * dfo_sources.ti_conc_02_ppm * dfo_sources.mass) / \
            np.sum(0.07437 * dfo_sources.ti_conc_02_ppm * dfo_sources.mass)
        f_03 = \
            (0.07437 * dfo_sources.ti_conc_03_ppm * dfo_sources.mass) / \
            np.sum(0.07437 * dfo_sources.ti_conc_03_ppm * dfo_sources.mass)

        e50ti_01[ii] = np.sum(f_01 * dfo_sources.e50ti_01)
        e50ti_02[ii] = np.sum(f_02 * dfo_sources.e50ti_02)
        e50ti_03[ii] = np.sum(f_03 * dfo_sources.e50ti_03)
        
    # Append to Dataframe
    dfo.loc[:,'ti_conc_01_ppm'] = ti_conc_01_ppm
    dfo.loc[:,'ti_conc_02_ppm'] = ti_conc_02_ppm
    dfo.loc[:,'ti_conc_03_ppm'] = ti_conc_03_ppm
    dfo.loc[:,'e50ti_01'] = e50ti_01
    dfo.loc[:,'e50ti_02'] = e50ti_02
    dfo.loc[:,'e50ti_03'] = e50ti_03
    
    # Return
    return dfo
