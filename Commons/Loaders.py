from Structs import Snapshot, Particle

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
    def __init__(self, nstep=1, ellipses=False):
        fname = "Outgasrun_%012d.dat" % nstep
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
                particle.id = float(line[1])                # -
                particle.m  = float(line[2])                # Msun
                particle.r  = float(line[3])                # 
                particle.x  = float(line[4])
                particle.y  = float(line[5])
                particle.z  = float(line[6])
                particle.vx = float(line[7])
                particle.vy = float(line[8])
                particle.vz = float(line[9])
                particle.cart2kep()
                if self.snapshot.ellipses:
                    particle.compute_ellipse()
                particles.append(particle)
            self.snapshot.particles = particles

class SSAscii(Loader):
    def __init__(self, nstep=1, ellipses=False):
        fname = "Out.%010d.dat" % nstep
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
                    self.snapshot.tout = float(line[0]) / 2. / pi       # TU -> yr
                    first = False
                particle = Particle()
                particle.id = float(line[1])
                particle.m  = float(line[2])
                particle.r  = float(line[3])
                particle.x  = float(line[4])
                particle.y  = float(line[5])
                particle.z  = float(line[6])
                particle.vx = float(line[7])
                particle.vy = float(line[8])
                particle.vz = float(line[9])
                particle.cart2kep()
                if self.snapshot.ellipses:
                    particle.compute_ellipse()
                particles.append(particle)
            self.snapshot.particles = particles

class SSBinary(Loader):
    pass
