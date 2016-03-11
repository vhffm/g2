"""
Compute Final Titanium Concentrations & e(50Ti/47Ti) Ratios.

Dirlist Format:
/zbox/data/volker/Debris/Runs/Chaos-EJS_Steep/gas_01/run_01
/zbox/data/volker/Debris/Runs/Chaos-EJS_Steep/gas_01/run_02
...
/zbox/data/volker/Debris/Runs/Chaos-CJS/gas_03/run_12
"""

import io_helpers as ioh
import isotope_helpers as iso
import sys
import pandas as pd


# List of Directories
if sys.stdin.isatty():
    print "!! No Directory List (Use Stdin)."
    sys.exit()
else:
    lines = sys.stdin.read().rstrip("\n").split("\n")
    dirs = []
    for line in lines:
        dirs.append(line)
    print "// Reading %i Directories" % len(dirs)

# Generate Tags & Basedirs
sim_tags = []
gas_tags = []
run_tags = []
basedirs = []

for line in lines:
    line_split = line.split('/')
    sim_tags.append(line_split[-3])
    gas_tags.append(line_split[-2])
    run_tags.append(line_split[-1])
    basedirs.append(line)

# Loop Simulations
for ibase, basedir in enumerate(basedirs):

    # Debug
    print "// Processing %s/%s/%s" % (sim_tags[ibase], \
                                      gas_tags[ibase], \
                                      run_tags[ibase])

    # Load Collisions
    fname = "%s/Collisions_%s.dat" % (basedir, run_tags[ibase])
    dfc = ioh.read_collisions_and_stack([fname])

    # Load IC Output
    fname = "%s/Out_%s_%012d.dat" % (basedir, run_tags[ibase], 0)
    dfo_t0 = ioh.read_output_and_stack([fname], frame='heliocentric')
    dfo_t0 = dfo_t0[dfo_t0.mass < 5.0]

    # Load Final Output
    fname = "%s/Out_%s_%012d.dat" % (basedir, run_tags[ibase], 9e9)
    dfo_tf = ioh.read_output_and_stack([fname], frame='heliocentric')
    dfo_tf = dfo_tf[dfo_tf.mass < 5.0]

    # Compute Concentrations
    dfo_tf = iso.compute_isotope_fractions(dfo_tf.copy(), \
                                           dfo_t0.copy(), \
                                           dfc.copy(), \
                                           showstep=False)

    # Add Sim/Gas/Run Identifiers to Dataframes
    dfo_tf['sim_tag'] = len(dfo_tf) * [ sim_tags[ibase] ]
    dfo_tf['gas_tag'] = len(dfo_tf) * [ gas_tags[ibase] ]
    dfo_tf['run_tag'] = len(dfo_tf) * [ run_tags[ibase] ]

    # Keep Interesting Columns
    dfo_tf = dfo_tf[['sim_tag', 'gas_tag', 'run_tag', \
                     'mass', 'a', \
                     'ti_conc_01_ppm', 'ti_conc_02_ppm', 'ti_conc_03_ppm', \
                     'e50ti_01', 'e50ti_02', 'e50ti_03']]

    # Concatenate to Output (Isotope) Dataframe
    if ibase == 0:
        dfo_iso = dfo_tf.copy()
    else:
        dfo_iso = pd.concat([dfo_iso, dfo_tf])

# Reset Indices
dfo_iso = dfo_iso.reset_index(drop=True)

# Outputting to CSV
print "// Writing to CSV"
dfo_iso.to_csv('titanium.csv', index=False)
