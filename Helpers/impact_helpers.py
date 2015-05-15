"""
Impact Simulation Helpers.
"""

import numpy as np
import constants as C


def calibrate_craters_n83(dfc_all, tt0_all, blk_all, target_pid=2):
    """
    Calibrate Collision Rate to Neukum (1983) Cratering Fits.
    Cf. Neukum+ 2001, Fig. 10 / Eqn. 5
    http://link.springer.com/article/10.1023/A:1011989004263

    @param: dfc_all - List of Collision Logs [List o/ Pandas Dataframe]
    @param: tt0_all - List of Initial Times (Years) [Numpy Float Array]
    @param: blk_all - List of Blacklisted Particles [List o/ Numpy Arrays]
    @param: target_pid - Particle ID of Earth (2 or 3) [Integer]
    @return: crc_all - Correction Factors [Numpy Float Array]
    """

    crc_all = []
    for idfc, dfc in enumerate(dfc_all):
        dfc = dfc[dfc.pidi==pid_target]

        # Time
        time = dfc[~dfc.pidj.isin(blk_all[idfc]) & \
                   (dfc.pidi==pid_target)].time[::-1] \
        time += t00_all[idfc] # Beginning of Sim
        time /= 1.0e9         # Yr => Gyr
        time += 0.115         # Gyr
        time += C.t0ss        # Beginning of Solar System
        time = np.asarray(time)

        # Cummulative Collision Ccount
        coll = np.cumsum(np.ones(len(dfc[~dfc.pidj.isin(blk_all[idfc]) & \
                                         (dfc.pidi==pid_target)])))
        coll /= C.Aearth

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
