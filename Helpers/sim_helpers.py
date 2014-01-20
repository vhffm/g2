"""
Simulation Helpers.
"""

import numpy as np

def mkparam(dt = None, \
            output_name = None, \
            energy_interval = 1, \
            coordinate_interval = 1, \
            integration_steps = None, \
            restart_step = None):
    """
    Generate Parameter File.
    """
    
    lines = []
    lines.append("Time step in days = %.2f\n" % dt)
    lines.append("Output name = %s\n" % output_name)
    lines.append("Energy output intervall = %i\n" % energy_interval)
    lines.append("Coordinates output intervall = %i\n" % coordinate_interval)
    lines.append("Number of outputs per intervall = 1\n")
    lines.append("Integration steps = %i\n" % integration_steps)
    lines.append("Central Mass = 1.0\n")
    lines.append("n1 = 3.0\n")
    lines.append("n2 = 0.4\n")
    lines.append("Input file = initial.dat\n")
    lines.append("Input file Format: << t i m r x y z vx vy vz Sx Sy Sz >>\n")
    lines.append("Default rho = 2.0\n")
    lines.append("Use Test Particles = 0\n")
    lines.append("Restart timestep = %i\n" % restart_step)
    lines.append("Minimum number of bodies = 1\n")
    lines.append("Order of integrator = 2\n")
    lines.append("aeGrid amin = 0\n")
    lines.append("aeGrid amax = 20\n")
    lines.append("aeGrid emin = 0\n")
    lines.append("aeGrid emax = 1\n")
    lines.append("aeGrid Na = 10\n")
    lines.append("aeGrid Ne = 10\n")
    lines.append("aeGrid Start Count = 10000000000000\n")
    lines.append("aeGrid name = A\n")
    lines.append("Gas dTau_diss = 3000000\n")
    lines.append("Gas alpha = 1\n")

    return lines
