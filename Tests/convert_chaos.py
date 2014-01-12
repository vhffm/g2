"""
Load XChaos. Read all vars. Bin lce/ltime. Write XChaos.
"""

from chaos_helpers import bin_lyapunov
import argparse
from time import gmtime, strftime
import numpy as np

# Parse Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--filename', default='XChaos.npz', \
                    help="Output File Name.")
args = parser.parse_args()

# Feedback
print "// Filename is %s" % args.filename

# Load Data
print "// (%s UTC) Loading Data" % strftime("%H:%M:%S", gmtime())
npz = np.load("%s" % args.filename)
lce = npz["lce"]; ds = npz["ds"]
istep0 = npz["istep0"]; nsteps = npz["nsteps"]
tout = npz["tout"]
c1pid0 = npz["c1pid0"]; c2pid0 = npz["c2pid0"]
c1a0 = npz["c1a0"]; c2a0 = npz["c2a0"]

# Bin Separations
print "// (%s UTC) Binning Separations" % strftime("%H:%M:%S", gmtime())
ds_mean, ds_median, ds_std, a0_bin_mids = bin_lyapunov(c1a0, ds)

# Bin Separation Ratios
print "// (%s UTC) Binning Separation Ratios" % strftime("%H:%M:%S", gmtime())
ds_ratio_mean, ds_ratio_median, ds_ratio_std, a0_bin_mids = \
    bin_lyapunov(c1a0, ds/ds[1,:])

# Bin Lyapunov Exponents
print "// (%s UTC) Binning LCEs" % strftime("%H:%M:%S", gmtime())
lce_mean, lce_median, lce_std, a0_bin_mids = bin_lyapunov(c1a0, lce)

# Bin Lyapunov Times
print "// (%s UTC) Binning Lyapunov Times" % strftime("%H:%M:%S", gmtime())
ltime_mean, ltime_median, ltime_std, a0_bin_mids = bin_lyapunov(c1a0, 1.0/lce)

# Save Relevant Arrays
print "// (%s UTC) Saving Data" % strftime("%H:%M:%S", gmtime())
np.savez("%s" % args.filename, \
    lce = lce, ds = ds, istep0 = istep0, nsteps = nsteps, tout = tout, \
    lce_mean = lce_mean, lce_median = lce_median, lce_std = lce_std, \
    ltime_mean = ltime_mean, ltime_median = ltime_median, \
    ltime_std = ltime_std, \
    ds_mean = ds_mean, ds_median = ds_median, ds_std = ds_std, \
    ds_ratio_mean = ds_ratio_mean, ds_ratio_median = ds_ratio_median, \
    ds_ratio_std = ds_ratio_std, \
    a0_bin_mids = a0_bin_mids, \
    c1pid0 = c1pid0, c2pid0 = c2pid0, c1a0 = c1a0, c2a0 = c2a0)

print "// Done"
