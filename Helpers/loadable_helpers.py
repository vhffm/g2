"""
Helpers to define what runs are loadable.
"""

def formation_sims():
    """
    Tracks what Formation runs are loadable (finished).
    Updated 06 Oct 2015.
    """
    
    # Set Directories
    basedir = "/zbox/data/volker/Debris/Runs"
    sims = [ "Chaos-VAN_ReRun/gas_00", \
             "Chaos-VAN_ReRun/gas_01", \
             "Chaos-VAN_ReRun/gas_03", \
             "Chaos-VAN_ReRun/gas_05", \
             "Chaos-VAN_Heavy/gas_01", \
             "Chaos-VAN_Heavy_4k/gas_01", \
             "Chaos-VAN_Steep/gas_01", \
             "Chaos-VAN_8192/gas_01_too_much_drag" ]

    # Valid Runs
    nruns = [ range(1,6+1),\
              range(1,6+1), range(1,6+1), range(1,6+1),\
              range(1,9+1), \
              range(1,5+1), \
              range(1,9+1),\
              range(1,12+1) ]

    # Return
    return basedir, sims, nruns


def chaos_sims():
    """
    Tracks what Chaos runs are loadable (finished).
    Updated 08 Oct 2015.
    """
    
    # Set Directories
    basedir = "/zbox/data/volker/Debris/Runs"
    sims = [ "Chaos-VAN/gas_01", \
             "Chaos-VAN_Steep/gas_01", \
             "Chaos-VAN_Heavy/gas_01", \
             "Chaos-EJS/gas_01", \
             "Chaos-CJS/gas_01", \
             "Chaos-EJS_Steep/gas_01", \
             "Chaos-EJS_Heavy/gas_01" ]

    # Valid Runs
    nruns = [ range(1,12+1),\
              range(1,12+1),\
              range(1,9+1),\
              range(1,12+1), \
              range(1,12+1), \
              range(1,6+1), \
              range(1,5+1) ]

    # Return
    return basedir, sims, nruns
