"""
Helpers to define what runs are loadable.
"""

def formation_sims():
    """
    Tracks what Formation runs are loadable (finished).
    Updated 28 Sep 2015
    """"
    
    # Set Directories
    basedir = "/zbox/data/volker/Debris/Runs"
    sims = [ "Chaos-VAN_ReRun/gas_00", \
             "Chaos-VAN_ReRun/gas_01", \
             "Chaos-VAN_ReRun/gas_03", \
             "Chaos-VAN_ReRun/gas_05", \
             "Chaos-VAN_Heavy/gas_01", \
             "Chaos-VAN_Heavy_4k/gas_01", \
             "Chaos-VAN_Steep/gas_01", \
             "Chaos-VAN_8192/gas_01" ]

    # Valid Runs
    nruns = [ range(1,6+1),\
              range(1,6+1), range(1,6+1), range(1,6+1),\
              range(1,9+1), \
              range(1,5+1), \
              range(1,9+1),\
              [ 1, 2, 3, 4, 5, 7, 8, 10, 11, 12 ] ]

    # Return
    return basedir, sims, nruns
