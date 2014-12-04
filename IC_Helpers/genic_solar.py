"""
Generate Initial Conditions for Solar System.
Based on Present-Day Solar System (Earth/Moon Barycenter).
Can Modify Outer Planets to EJS/CJS Orbits.
Output Format Is Default Genga IC.
"""

import numpy as np
import argparse
import kepler_helpers as kh

# Parse arguments
parser = argparse.ArgumentParser()
group1 = parser.add_mutually_exclusive_group()
group1.add_argument('--ejs', action='store_true', help="EJS Configuration.")
group1.add_argument('--cjs', action='store_true', help="CJS Configuration.")
group1.add_argument('--nice1', action='store_true', help="Nice 1 Model.")
group1.add_argument('--nice2', action='store_true', help="Nice 2 Model.")
group1.add_argument('--nice3', action='store_true', help="Nice 3 Model.")
group1.add_argument('--nice4', action='store_true', help="Nice 4 Model.")
group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument('--all', action='store_true', help="All Planets.")
group2.add_argument('--outer', action='store_true', help="Outer Planets.")
group2.add_argument('--inner', action='store_true', help="Inner Planets.")
group2.add_argument('--js', action='store_true', help="Jupiter/Saturn Only.")
args = parser.parse_args()

# NASA Horizon Query - 01-01-2014
# Earth Without Moon
# Deprecated
# mercury = "0 0 1.6515006786989092e-07 1.6310392545626536e-05 +1.1972692892595370e-01 -4.3426096069639630e-01 -4.6466929699041107e-02 +2.1482009903780699e-02 +8.9119212513564267e-03 -1.2427988748894680e-03 0 0 0"
# venus = "0 1 2.4468352521240763e-06 4.0455121182840900e-05 -4.9881994854843653e-02 +7.1760434138711926e-01 +1.2711854910906791e-02 -2.0246486048184902e-02 -1.5114700566799760e-03 +1.1477577544001651e-03 0 0 0"
# earth = "0 2 3.0023628776833749e-06 4.2587504470568303e-05 -1.7558926311807671e-01 +9.6752446474730014e-01 -2.9820874849228900e-05 -1.7207909676081830e-02 -3.1365065792022962e-03 +1.1244539739591660e-07 0 0 0"
# mars = "0 3 3.2125081695239055e-07 2.2660750299046701e-05 -1.5124387605410869e+00 +6.9681542097585980e-01 +5.1724065828569588e-02 -5.3309348873500006e-03 -1.1514537719341871e-02 -1.1040581507850680e-04 0 0 0"
# jupiter = "0 4 9.5420039213714751e-04 4.6732616936774455e-04 -1.3307968217788140e+00 +5.0187260363223007e+00 +8.9354534335869748e-03 -7.3915532677911646e-03 -1.5771896615651140e-03 +1.7200784821425221e-04 0 0 0"
# saturn = "0 5 2.8570710371524815e-04 3.8925687652332965e-04 -6.8851925827531151e+00 -7.0754782545976083e+00 +3.9713562554594872e-01 +3.6895716420435351e-03 -3.9070640598097067e-03 -7.8752021502553231e-05 0 0 0"
# uranus = "0 6 4.3642853551857632e-05 1.6953449825499186e-04 +1.9645011082335920e+01 +3.9225997140108890e+00 -2.3986827787041001e-01 -8.0447801422034405e-04 +3.6714036954299551e-03 +2.4141556323116709e-05 0 0 0"
# neptune = "0 7 5.1480569101603742e-05 1.6458790379443301e-04 +2.7063448807961450e+01 -1.2891858367498200e+01 -3.5810433631098909e-01 +1.3227724210187141e-03 +2.8500466203526960e-03 -8.9055141344433750e-05 0 0 0"
# pluto = "0 8 6.5808657181639943e-09 7.7608056333903304e-06 +6.2567672814280284e+00 -3.1924879744181990e+01 +1.6063792641261310e+00 +3.1299238024423749e-03 -3.0015467471074471e-05 -8.9263215128274780e-04 0 0 0"

# Joachim's ICs - ??-??-????
# Earth = Earth/Moon Barycenter
mercury = "0 0 1.6515006786989092e-07 1.6310392545626536e-05  3.704735169720974E-02 -4.529211095852149E-01 -4.090255306376755E-02 1.301692241729103 0.21718328966069625 -0.10173333623267612 0 0 0"
venus = "0 1 2.4468352521240763e-06 4.04551211828409e-05 4.272157290820016E-01 -5.835752726996720E-01 -3.279422047835795E-02 0.9430994394414915 0.6869103343210567 -0.045042426747977045 0 0 0"
earth = "0 2 3.002362877683375e-06 4.25875044705683e-05 -9.948486972722731E-01 4.564231864395614E-02 -6.099525188647536E-05 -0.05755930245755808 -1.0030463160527843 2.5265763469203405e-05 0 0 0"
mars = "0 3 3.2125081695239055e-07 2.26607502990467e-05 -1.093539305796724E+00 1.240381444357973E+00 5.266905384900308E-02 -0.5789102231935592 -0.4698447773869488 0.004371971509388125 0 0 0"
jupiter = "0 4 0.0009542003921371475 0.00046732616936774455 7.199075962861715E-01 -5.164765414047316E+00 5.281301305052329E-03 0.4290585620283851 0.08134728734237137 -0.009940783595198496 0 0 0"
saturn = "0 5 0.00028570710371524815 0.00038925687652332965 -8.469664737705321E+00 3.804527121928150E+00 2.708474727487031E-01 -0.15016129989332302 -0.2965903629571777 0.011135663432612576 0 0 0"
uranus = "0 6 4.364285355185763e-05 0.00016953449825499186 1.970001443062262E+01 -3.956376098536538E+00 -2.699288868040702E-01 0.04335170238397214 0.2135087007559093 0.00023187752934430736 0 0 0"
neptune = "0 7 5.148056910160374e-05 0.000164587903794433 2.361441531200179E+01 -1.856288724958460E+01 -1.619425696998957E-01 0.1115723331088743 0.1445374992111762 -0.005547767271146534 0 0 0"
pluto = "0 8 6.580865718163994e-09 7.76080563339033e-06 -4.656585770964581E-01 -3.123136435608064E+01 3.476634650132539E+00 0.18578751327844598 -0.03702089125652011 -0.04977940553539111 0 0 0"

# Adjust Orbital Parameters
if args.cjs or args.ejs or args.nice1 or args.nice2:
    # CJS, EJS Models
    # Cf. Raymond 2006, Morishima 2010
    if args.cjs:
        anew = np.array([5.45, 8.18])
        enew = np.array([0.00, 0.00])
        inew = np.array([0.00, 0.00])
    if args.ejs:
        anew = np.array([5.200, 9.550])
        enew = np.array([0.048, 0.056])
        inew = np.array([0.020, 0.000])

    # Nice1 Model
    # Cf. Tsiganis 2005 (Text, Figure 1)
    # The Nice model ICs have different masses for Ice Giants. Let's ignore.
    if args.nice1:
        # Jupiter, Saturn, Uranus, Neptune
        anew = np.array([5.45, 8.65, 12.00, 15.25])
        enew = np.array([0.001, 0.001, 0.001, 0.001])
        inew = np.array([0.001, 0.001, 0.001, 0.001])

    # Nice2 Model
    # Cf. Levison 2011+ (Table 1), Morbidelli+ 2007
    # Cf. http://www.lpi.usra.edu/meetings/lpsc2013/eposter/2772.pdf
    # The Nice model ICs have different masses for Ice Giants. Let's ignore.
    if args.nice2:
        anew = np.array([5.42, 7.32, 9.61, 11.67])
        enew = np.array([0.0044, 0.017, 0.053, 0.011])
        inew = np.array([0.016, 0.016, 0.044, 0.029])

    # Nice3 Model
    # Cf. Figure 1 Caption in Gomes+ 2005
    if args.nice3:
        anew = np.array([5.45, 8.18, 11.5, 14.2])
        enew = np.array([0.001, 0.001, 0.001, 0.001])
        inew = np.array([0.001, 0.001, 0.001, 0.001])

    # Nice4 Model
    # Like Nice3, but Jupiter and Saturn P1/P2 ~ 1.98
    if args.nice4:
        anew = np.array([5.45, 8.6, 11.5, 14.2])
        enew = np.array([0.0001, 0.0001, 0.0001, 0.0001])
        inew = np.array([0.0001, 0.0001, 0.0001, 0.0001])

    # What Planets To Fix?
    if args.ejs or args.cjs:
        plist = [ jupiter, saturn ]
    elif args.nice1 or args.nice2 or args.nice3 or args.nice4:
        plist = [ jupiter, saturn, uranus, neptune ]

    for iplanet, planet in enumerate(plist):
        line = planet.strip().split()
        iloc = int(line[1])
        mloc = float(line[2])
        rloc = float(line[3])
        xold = np.array([ float(line[4]), float(line[5]), float(line[6]) ])
        vold = np.array([ float(line[7]), float(line[8]), float(line[9]) ])
        _, _, _, Omega, omega, M0 = kh.cart2kep(xold, vold, mloc)
        xnew, vnew = kh.kep2cart(anew[iplanet], enew[iplanet], inew[iplanet], \
                                 Omega, omega, M0, mloc)
        # <<< t i m r x y z vx vy vz Sx Sy Sz >>>
        line = "0 %i %.19f %.19f %.19f %.19f %.19f %.19f %.19f %.19f 0.0 0.0 0.0" % (iloc, mloc, rloc, xnew[0], xnew[1], xnew[2], vnew[0], vnew[1], vnew[2])
        if iplanet == 0:
            jupiter = line
        elif iplanet == 1:
            saturn = line
        elif iplanet == 2:
            uranus = line
        elif iplanet == 3:
            neptune = line

# Cat
if args.all:
    ss = [ mercury, venus, earth, mars, \
           jupiter, saturn, uranus, neptune, pluto ]
if args.js:
    ss = [ jupiter, saturn ]
if args.outer:
    ss = [ jupiter, saturn, uranus, neptune, pluto ]
if args.inner:
    ss = [ mercury, venus, earth, mars ]

# Stdout
for planet in ss:
    print planet
