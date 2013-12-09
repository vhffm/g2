"""
Chop up particle files into blocks of args.blocklen.
Particle file must be in the form < t i m r x y z vx vy vz Sx Sy Sz >.
"""

import numpy as np
from time import gmtime, strftime
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--filename', type=str, \
                    help="Input Filename.", required=True)
parser.add_argument('--blocklen', type=int, default=100000, \
                    help="Block Length.")
args = parser.parse_args()

# Read test particles
print "// (%s UTC) Reading Particle File (%s)" % \
    (strftime("%H:%M:%S", gmtime()), args.filename)
with open(args.filename, 'r') as f:
    lines = f.readlines()
    nlines = len(lines)
    print "// (%s UTC) %i Particles Found" % \
        (strftime("%H:%M:%S", gmtime()), nlines)

    t = np.zeros(nlines)
    pid = np.zeros(nlines, dtype='int')
    m = np.zeros(nlines)
    r = np.zeros(nlines)
    x = np.zeros(nlines)
    y = np.zeros(nlines)
    z = np.zeros(nlines)
    vx = np.zeros(nlines)
    vy = np.zeros(nlines)
    vz = np.zeros(nlines)
    Sx = np.zeros(nlines)
    Sy = np.zeros(nlines)
    Sz = np.zeros(nlines)

    for iline, line in enumerate(lines):
        line = line.strip().split()
        t[iline] = float(line[0])
        pid[iline] = int(line[1])
        m[iline] = float(line[2])
        r[iline] = float(line[3])
        x[iline] = float(line[4])
        y[iline] = float(line[5])
        z[iline] = float(line[6])
        vx[iline] = float(line[7])
        vy[iline] = float(line[8])
        vz[iline] = float(line[9])
        Sx[iline] = float(line[10])
        Sy[iline] = float(line[11])
        Sz[iline] = float(line[12])

    # Compute (squared) current orbital distance
    d2 = np.sqrt(x**2. + y**2. + z**2.)

    # Indices that sort the arrays
    iisorted = np.argsort(d2)

    # Sort arrays
    t = t[iisorted]
    pid = pid[iisorted]
    m = m[iisorted]
    r = r[iisorted]
    x = x[iisorted]
    y = y[iisorted]
    z = z[iisorted]
    vx = vx[iisorted]
    vy = vy[iisorted]
    vz = vz[iisorted]
    Sx = Sx[iisorted]
    Sy = Sy[iisorted]
    Sz = Sz[iisorted]

# Write output files
iline = 0; iblock = 1; blocklen = args.blocklen
print "// (%s UTC) Writing Chopped Particle Files" % strftime("%H:%M:%S", gmtime())
print "//                Blocks of %i Particles" % blocklen
while iline < nlines:
    fout = "testparticles_%05d.dat" % iblock
    print "// (%s UTC) Writing %s" % (strftime("%H:%M:%S", gmtime()), fout)
    with open(fout, 'w') as f:
        for jline in range(iline, iline + blocklen):
            # Prevent iline + blocklen > nlines
            if jline == nlines:
                break
            linestr  = "%.16e " % t[jline]
            linestr += "%07d " % pid[jline]
            linestr += "%.8e " % m[jline]
            linestr += "%.8e " % r[jline]
            linestr += "%+.8e " % x[jline]
            linestr += "%+.8e " % y[jline]
            linestr += "%+.8e " % z[jline]
            linestr += "%+.8e " % vx[jline]
            linestr += "%+.8e " % vy[jline]
            linestr += "%+.8e " % vz[jline]
            linestr += "%+.8e " % Sx[jline]
            linestr += "%+.8e " % Sy[jline]
            linestr += "%+.8e " % Sz[jline]
            linestr += "\n"
            f.write(linestr)
            iline += 1
        iblock += 1
print "// (%s UTC) Done" % strftime("%H:%M:%S", gmtime())
