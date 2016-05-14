"""
Generate Initial Conditions for Solar System.
Based on Present-Day Solar System (Earth/Moon Barycenter).
Two Options For Present-Day SS. Use NASA Horizons Output When In Doubt.
Can Modify Outer Planets to EJS/CJS Orbits.
Output Format Is Default Genga IC.
"""

import numpy as np
import argparse
import kepler_helpers as kh
import constants as C
import ic_helpers as ih

# Parse arguments
parser = argparse.ArgumentParser()
group1 = parser.add_mutually_exclusive_group()
group1.add_argument('--ejs', action='store_true', help="EJS Configuration.")
group1.add_argument('--cjs', action='store_true', help="CJS Configuration.")
group1.add_argument('--nice1', action='store_true', help="Nice 1 Model.")
group1.add_argument('--nice2', action='store_true', help="Nice 2 Model.")
group1.add_argument('--nice3', action='store_true', help="Nice 3 Model.")
group1.add_argument('--nice4', action='store_true', help="Nice 4 Model.")
group1.add_argument('--morby', action='store_true', \
                    help="Alessandro Morbidelli Mail.")
group2 = parser.add_mutually_exclusive_group(required=True)
group2.add_argument('--all', action='store_true', help="All Planets.")
group2.add_argument('--outer', action='store_true', help="Outer Planets.")
group2.add_argument('--inner', action='store_true', help="Inner Planets.")
group2.add_argument('--giants', action='store_true', \
                    help="Giant (Gas/Ice) Planets.")
group2.add_argument('--js', action='store_true', help="Jupiter/Saturn Only.")
group3 = parser.add_mutually_exclusive_group(required=True)
group3.add_argument('--solar1', action='store_true', \
                    help="Base Solar System: Joachim")
group3.add_argument('--solar2', action='store_true', \
                    help="Base Solar System: NASA Horizron Query 2014-01-01.")
group3.add_argument('--solar3', action='store_true', \
                    help="Base Solar System: NASA Horizron Query 2015-03-19.")
args = parser.parse_args()

# Joachim's ICs
# Origin = ???
# Epoch = ???
# Earth = Earth/Moon Barycenter
if args.solar1:
    mercury = "0 1 1.6515006786989092e-07 1.6310392545626536e-05  3.704735169720974E-02 -4.529211095852149E-01 -4.090255306376755E-02 1.301692241729103 0.21718328966069625 -0.10173333623267612 0 0 0"
    venus = "0 2 2.4468352521240763e-06 4.04551211828409e-05 4.272157290820016E-01 -5.835752726996720E-01 -3.279422047835795E-02 0.9430994394414915 0.6869103343210567 -0.045042426747977045 0 0 0"
    earth = "0 3 3.002362877683375e-06 4.25875044705683e-05 -9.948486972722731E-01 4.564231864395614E-02 -6.099525188647536E-05 -0.05755930245755808 -1.0030463160527843 2.5265763469203405e-05 0 0 0"
    mars = "0 4 3.2125081695239055e-07 2.26607502990467e-05 -1.093539305796724E+00 1.240381444357973E+00 5.266905384900308E-02 -0.5789102231935592 -0.4698447773869488 0.004371971509388125 0 0 0"
    jupiter = "0 5 0.0009542003921371475 0.00046732616936774455 7.199075962861715E-01 -5.164765414047316E+00 5.281301305052329E-03 0.4290585620283851 0.08134728734237137 -0.009940783595198496 0 0 0"
    saturn = "0 6 0.00028570710371524815 0.00038925687652332965 -8.469664737705321E+00 3.804527121928150E+00 2.708474727487031E-01 -0.15016129989332302 -0.2965903629571777 0.011135663432612576 0 0 0"
    uranus = "0 7 4.364285355185763e-05 0.00016953449825499186 1.970001443062262E+01 -3.956376098536538E+00 -2.699288868040702E-01 0.04335170238397214 0.2135087007559093 0.00023187752934430736 0 0 0"
    neptune = "0 8 5.148056910160374e-05 0.000164587903794433 2.361441531200179E+01 -1.856288724958460E+01 -1.619425696998957E-01 0.1115723331088743 0.1445374992111762 -0.005547767271146534 0 0 0"
    pluto = "0 9 6.580865718163994e-09 7.76080563339033e-06 -4.656585770964581E-01 -3.123136435608064E+01 3.476634650132539E+00 0.18578751327844598 -0.03702089125652011 -0.04977940553539111 0 0 0"

# NASA Horizon Query
# Origin = Sun, Body Center
# Epoch = 01-01-2014, 00:00 (Coordinate Time)
# Earth = Earth (No Moon)
if args.solar2:
    plist, _ = ih.Solar2(epoch='2014-01-01')
    # Spaghetti code ahead
    for iplanet, planet in enumerate(plist):
        if iplanet == 0:
            mercury = planet
        elif iplanet == 1:
            venus = planet
        elif iplanet == 2:
            earth = planet
        elif iplanet == 3:
            mars = planet
        elif iplanet == 4:
            jupiter = planet
        elif iplanet == 5:
            saturn = planet
        elif iplanet == 6:
            uranus = planet
        elif iplanet == 7:
            neptune = planet
        elif iplanet == 8:
            pluto = planet

# NASA Horizon Query
# Origin = Sun, Body Center
# Epoch = 19-03-2015, 00:00 (Barycentric Dynamical Time)
# Earth = Earth (No Moon)
if args.solar3:
    plist, _ = ih.Solar2(epoch='2015-03-19')
    # Spaghetti code ahead
    for iplanet, planet in enumerate(plist):
        if iplanet == 0:
            mercury = planet
        elif iplanet == 1:
            venus = planet
        elif iplanet == 2:
            earth = planet
        elif iplanet == 3:
            mars = planet
        elif iplanet == 4:
            jupiter = planet
        elif iplanet == 5:
            saturn = planet
        elif iplanet == 6:
            uranus = planet
        elif iplanet == 7:
            neptune = planet
        elif iplanet == 8:
            pluto = planet

# Adjust Orbital Parameters
if args.cjs or args.ejs or \
    args.nice1 or args.nice2 or args.nice3 or args.nice4 or \
    args.morby:

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

    # ICs for 1 Gyr Stable Run in Morbidelli+ 2007 / Levison+ 2011
    # Reported Stable for 4 Gyr in Tests
    # Source: E-Mail
    # Date: Tue, 17 Feb 2015 11:28:50 +0100
    # From: alessandro morbidelli <morby@oca.eu>
    # To: Volker Hoffmann <volker@physik.uzh.ch>, hal@boulder.swri.edu
    if args.morby:
        # semi major axis (AU), eccentricity, inclination (rad), longitude of the node (rad), argument of pericenter (rad), mean anomaly (rad), mass (n units where the Sun's mass is (2PI)^2.
        #   0.5429655E+01  0.4890734E-02  0.6005074E-03  0.2099288E+01  0.6161135E+01  0.1208643E+00  0.3947842E-01
        #   0.7298669E+01  0.9955255E-02  0.1175241E-03  0.2056477E+01  0.2942986E+01  0.3344157E+01  0.1121187E-01
        #   0.9640065E+01  0.4988758E-01  0.9948838E-03  0.1854808E+01  0.2748019E+01  0.3574791E+01  0.1776529E-02
        #   0.1161223E+02  0.1060971E-01  0.1614985E-03  0.2498409E+01  0.4613698E+01  0.4832179E+01  0.1776529E-02
        anew = np.array([0.5429655E+01, 0.7298669E+01, \
                         0.9640065E+01, 0.1161223E+02])
        enew = np.array([0.4890734E-02, 0.9955255E-02, \
                         0.4988758E-01, 0.1060971E-01])
        inew = np.array([0.6005074E-03, 0.1175241E-03, \
                         0.9948838E-03, 0.1614985E-03])
        Omega = np.array([0.2099288E+01, 0.2056477E+01, \
                          0.1854808E+01, 0.2498409E+01])
        omega = np.array([0.6161135E+01, 0.2942986E+01, \
                          0.2748019E+01, 0.4613698E+01])
        M0 = np.array([0.1208643E+00, 0.3344157E+01, \
                          0.3574791E+01, 0.4832179E+01])
        mnew = np.array([0.3947842E-01, 0.1121187E-01, \
                         0.1776529E-02, 0.1776529E-02])
        mnew /= C.twopi**2.0

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
        enew = np.array([0.0001, 0.0001, 0.0001, 0.0001])
        inew = np.array([0.0001, 0.0001, 0.0001, 0.0001])

    # Nice4 Model
    # Like Nice3, but Jupiter and Saturn P1/P2 ~ 1.98
    if args.nice4:
        anew = np.array([5.45, 8.6, 11.5, 14.2])
        enew = np.array([0.0001, 0.0001, 0.0001, 0.0001])
        inew = np.array([0.0001, 0.0001, 0.0001, 0.0001])

    # What Planets To Fix?
    if args.ejs or args.cjs:
        plist = [ jupiter, saturn ]
    elif args.nice1 or args.nice2 or args.nice3 or args.nice4 or args.morby:
        plist = [ jupiter, saturn, uranus, neptune ]

    for iplanet, planet in enumerate(plist):
        line = planet.strip().split()
        iloc = int(line[1])
        mloc = float(line[2])
        rloc = float(line[3])
        # Nice{1,2,3,4} & EJS => Keep Angles
        if not args.morby:
            xold = np.array([ float(line[4]), float(line[5]), float(line[6]) ])
            vold = np.array([ float(line[7]), float(line[8]), float(line[9]) ])
            _, _, _, Omega, omega, M0 = kh.cart2kep(xold, vold, mloc)
            xnew, vnew = \
                kh.kep2cart(anew[iplanet], enew[iplanet], inew[iplanet], \
                            Omega, omega, M0, mloc)
        # Morby ICs
        else:
            xnew, vnew = kh.kep2cart(anew[iplanet], \
                                     enew[iplanet], \
                                     inew[iplanet], \
                                     Omega[iplanet], \
                                     omega[iplanet], \
                                     M0[iplanet], \
                                     mnew[iplanet])
            mloc = mnew[iplanet]
        # <<< t i m r x y z vx vy vz Sx Sy Sz >>>
        line = "0.0 %05d %.16e %.16e " % (iloc, mloc, rloc)
        line += "%+.16e %+.16e %+.16e " % (xnew[0], xnew[1], xnew[2])
        line += "%+.16e %+.16e %+.16e " % (vnew[0], vnew[1], vnew[2])
        line += "0.0 0.0 0.0"
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
if args.giants:
    ss = [ jupiter, saturn, uranus, neptune ]
if args.inner:
    ss = [ mercury, venus, earth, mars ]

# Stdout
for planet in ss:
    print planet
