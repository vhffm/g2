"""
Loading Routines.

NB: For G=M=1, 1 Yr = 2 Pi.
    Hence the velocity (1 Day = 2 pi / 365.25 ~ 0.0172020989)
"""

from Structs import Snapshot, Particle
import numpy as np

class Loader():
    def __init__(self, nstep, fname, ellipses):
        self.nstep = nstep
        self.fname = fname
        self.snapshot = Snapshot()
        self.snapshot.nstep = nstep
        self.snapshot.ellipses = ellipses

class GengaIC(Loader):
    pass

class GengaOut(Loader):
    def __init__(self, nstep=1, ellipses=False, run_name="gasrun"):
        fname = "Out%s_%012d.dat" % (run_name, nstep)
        Loader.__init__(self, nstep, fname, ellipses)

    def load(self):
        with open(self.fname, 'r') as f:
            lines = f.readlines()
            first = True
            particles = []
            self.snapshot.nparticles = len(lines)
            for line in lines:
                line = line.strip().split(" ")
                if first:
                    self.snapshot.tout = float(line[0])     # yr
                    first = False
                particle = Particle()
                particle.id = float(line[1])    # -
                particle.m  = float(line[2])    # Msun
                particle.r  = float(line[3])    # AU
                particle.x  = float(line[4])    # AU
                particle.y  = float(line[5])    # AU
                particle.z  = float(line[6])    # AU
                particle.vx = float(line[7])    # AU/day/0.0172020989 = AU/yr
                particle.vy = float(line[8])    # AU/day/0.0172020989 = AU/yr
                particle.vz = float(line[9])    # AU/day/0.0172020989 = AU/yr
                particle.cart2kep()
                if self.snapshot.ellipses:
                    particle.compute_ellipse()
                particles.append(particle)
            self.snapshot.particles = particles

class SSAscii(Loader):
    def __init__(self, nstep=1, ellipses=False):
        fname = "Out.%012d.dat" % nstep
        Loader.__init__(self, nstep, fname, ellipses)

    def load(self):
        with open(self.fname, 'r') as f:
            lines = f.readlines()
            first = True
            particles = []
            self.snapshot.nparticles = len(lines)
            for line in lines:
                line = line.strip().split(" ")
                if first:
                    self.snapshot.tout = float(line[0]) / 2. / np.pi    # TU -> yr
                    first = False
                particle = Particle()
                particle.id = float(line[1])    # -
                particle.m  = float(line[2])    # Msun
                particle.r  = float(line[3])    # AU
                particle.x  = float(line[4])    # AU
                particle.y  = float(line[5])    # AU
                particle.z  = float(line[6])    # AU
                particle.vx = float(line[7])    # AU/day/0.0172020989 = AU/yr
                particle.vy = float(line[8])    # AU/day/0.0172020989 = AU/yr
                particle.vz = float(line[9])    # AU/day/0.0172020989 = AU/yr
                particle.cart2kep()
                if self.snapshot.ellipses:
                    particle.compute_ellipse()
                particles.append(particle)
            self.snapshot.particles = particles

class SSBinary(Loader):
    pass
