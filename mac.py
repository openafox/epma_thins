#!/usr/bin/env python
"""-------------------------------MAC----------------------------------

a) calculation of the parameters used to calculate the
   mass absortion coefficients
b) calculation of mass absorption coefficients
   heinrich, proceedings icxom-11 1986.
10/88 r.a. waldo

Notes:
- will edit and check once I get access to the paper
- likley needs to be updated to newer MAC Caculations ?? research needed

Keyword arguments:
z -- Atomic Number
mass -- Atomic Mass
xray -- Element X-Ray emmision
Return:
xmu -- ??Mass Absorption Coefficient??
iflag -- ??
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import numpy as np
from basictools import get_data


def mac(z, mass, xray):
    """ This calc needs to be checked.  See above refs"""
    iflag = 1
    cutoff = 0
    e = 0
    # calculate 'c' parameter
    if xray > get_data(z, 'K'):
        ie = 1
        if z < 6:
            c = -2.87536e-4 + 1.808599e-3 * z
        else:
            c = (5.253e-3 +
                 1.33257e-3 * z +
                 -7.5937e-5 * z**2 +
                 1.69357e-6 * z**3 +
                 -1.3975e-8 * z**4)
    elif xray < get_data(z, 'L3'):
        c = (-9.24e-5 +
             1.41478e-4 * z +
             -5.24999e-6 * z**2 +
             9.85296e-8 * z**3 +
             -9.07306e-10 * z**4 +
             3.19254e-12 * z**5)
        if xray > get_data(z, 'L1'):
            ie = 2
            c = c
        elif xray > get_data(z, 'L2'):
            ie = 3
            c = c * 0.858
        else:
            ie = 4
            c = c*(0.8933 - z*8.29e-3 + z**2*6.38e-5)
    elif xray > get_data(z, 'M1'):
        ie = 5
        if z < 30:
            c = (1.889757e-2 +
                 -1.8517159e-3 * z +
                 6.9602789e-5 * z**2 +
                 -1.1641145e-6 * z**3 +
                 7.2773258e-9 * z**4)
        else:
            c = (3.0039e-3 +
                 -1.73663566e-4 * z +
                 4.0424792e-6 * z**2 +
                 -4.0585911e-8 * z**3 +
                 1.497763e-10 * z**4)
    elif xray > get_data(z, 'M5'):
        c1 = (7.7708e-5 +
              -7.83544e-6 * z +
              2.209365e-7 * z**2 +
              -1.29086e-9 * z**3)
        c2 = (1.406 +
              0.0162 * z +
              -6.561e-4 * z**2 +
              4.865e-6 * z**3)
        c3 = (0.584 +
              0.01955 * z +
              -1.285e-4 * z**2)
        c4 = (1.082 +
              1.366e-3 * z)
        c5 = (1.6442 +
              -0.0480 * z +
              4.0664e-4 * z**2)
        if xray > get_data(z, 'M2'):
            ie = 6
            c = c1*c2*c3
        elif xray > get_data(z, 'M3'):
            ie = 7
            c = c1*c2*c4
        elif xray > get_data(z, 'M4'):
            ie = 8
            c = c1*c2*0.95
        else:
            ie = 9
            c = c1*c2*c5
    else:
        ie = 10
        c = 1.08*(4.3156e-3 +
                  -1.4653e-4 * z +
                  1.707073e-6 * z**2 +
                  -6.69827e-9 * z**3)
        if xray < get_data(z, 'N1'):
            cutoff = ((0.252*z - 31.1812)*z + 1042.)/1000.0
            e = cutoff
            ie = 11

#      calculate 'n' parameter
    if xray > get_data(z, 'K'):
        if z < 6:
            n = (3.34745 +
                 0.02652873 * z +
                 -0.01273815 * z**2)
        else:
            n = (3.112 +
                 -0.0121 * z)
    elif xray > get_data(z, 'L3'):
        n = (2.7575 +
             1.889e-3 * z +
             -4.982e-5 * z**2)
    elif xray > get_data(z, 'M1'):
        n = (0.5385 +
             0.084597 * z +
             -1.08246e-3 * z**2 +
             4.4509e-6 * z**3)
    elif xray > get_data(z, 'M5'):
        n = 3.0 - 0.004*z
    else:
        n = 0.3736 + 0.02401*z

#      calculate 'a' parameter
    if xray > get_data(z, 'K'):
        if z < 6:
            a = (24.4545 +
                 155.6055 * z +
                 -14.15422 * z**2)
        else:
            a = (47.0 * z +
                 6.52 * z**2 +
                 -0.152624 * z**3)
    elif xray > get_data(z, 'L3'):
        a = (17.8096 * z +
             0.067429 * z**2 +
             0.01253775 * z**3 +
             -1.16286e-4 * z**4)
    elif xray > get_data(z, 'M1'):
        a = (10.2575657 * z +
             -0.822863477 * z**2 +
             2.63199611e-2 * z**3 +
             -1.8641019e-4 * z**4)
    elif xray > get_data(z, 'M5'):
        a = (4.62 * z +
             -0.04 * z**2)
    else:
        a = (19.64 * z +
             -0.61239 * z**2 +
             5.39309e-3 * z**3)

#      calculate 'b' parameter
    if xray > get_data(z, 'K'):
        if z < 6:
            b = -103. + 18.2*z
        else:
            b = 0.
    elif xray > get_data(z, 'L3'):
        b = 0.
    elif xray > get_data(z, 'M1'):
        if z < 61:
            b = (5.654 * z +
                 -0.536839169 * z**2 +
                 0.018972278 * z**3 +
                 -1.683474e-4 * z**4)
        else:
            b = (-1232.4022 * z +
                 51.114164 * z**2 +
                 -0.699473097 * z**3 +
                 3.1779619e-3 * z**4)
    elif xray > get_data(z, 'M5'):
        b = (2.51 +
             (-0.052) * z +
             3.78e-4 * z**2) * get_data(z, 'M4') * 1000.
    else:
        b = -113. + 4.5*z

    if xray > get_data(z, 'N1'):
        qq = (-xray*1000.+b)/a
        if qq > 88.:
            qq = 88.
        elif qq < -88.:
            qq = -88.

        xmu = c * z**4 / mass * ((12.397 / xray)**n) * (1. - np.exp(qq))
        if xmu < 0.0:
            iflag = 8
    else:
        xmu = (((12.397 / xray)**n) * c * z**4 / mass * (xray - cutoff) /
               (1.08 * get_data(z, 'N1')))
    if xray < 1.1 * cutoff:
        iflag = 3
    if xray - e < 0.02 and xray - e > -0.005:
        iflag = 2
    if ie == 9:
        iflag = 5
    if ie == 11:
        iflag = 4
    if ie == 11 and (xray - e < 0.02 and xray - e > -0.005):
        iflag = 6
    if ie == 9 and (xray - e < 0.02 and xray - e > -0.005):
        iflag = 7
    if xmu < 0.0 and iflag == 4:
        iflag = 9

    return (xmu, iflag)

if __name__ == '__main__':
    print mac(12, 24.3, 1.25)
