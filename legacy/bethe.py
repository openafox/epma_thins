#!/usr/bin/env python
"""This is my doc string.

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import numpy as np


def decel_pap(samp, el_p, line=None):
    """Returns Total Trajectory
        or Deceleration Factor if line is provided
        el_p -- primary element ie element of intrest of excited element
        [1] pages 34-6
    """
    # Overvoltage Ratio
    el_i=el_p
    u0 = el_p.volt/(el_p.xray/1000.)
    mavg = (el_i.z * el_i.c1) / el_i.mass
    j_i = el_i.z*(10.04 + 8.25*np.exp(-1.0*el_i.z/11.22))/1000.  # keV

    # [1]p35e8
    p_k = [0.78,
           0.1,
           (-1.*(0.5 - 0.25*j_i))]

    d_k = [6.6e-6,
           1.12e-5*(1.35 - 0.45*(j_i**2)),
           2.2e-6/j_i]

    con = ((1.602e-19)**2*(6.023e23)/(8*np.pi*(8.854e-12)**2))/100
    if line == 'dE/dps':
        # find just the average energy loss (for checking) [1]p34e5
        pap_dec = 0
        f_V = 0
        for k in range(0, 3):
            print(d_k[k])
            print(p_k[k])
            f_V += (d_k[k]*(el_p.volt/j_i)**(p_k[k]))
        pap_dec = (mavg/j_i)*(1.0/(f_V))

        # Bethe
        con = ((1.602e-19)**2*(6.023e23)/(8*np.pi*(8.854e-12)**2))/100
        # Con ~ 78500 -

        bethe_dec2 =(con/(el_p.volt)*
                    mavg*np.log(1.166*el_p.volt/j_i))
        return pap_dec, con, bethe_dec2


if __name__ == '__main__':
    # Need to set this up to test the pap - getting close to test 2/22/16
    from analysis_sample import AnalysisSample
    from atomic_element import AtomicElement as AtEl
    from film_layer import FilmLayer as FL
    import matplotlib.pyplot as plt
    el_p = AtEl('Cu', 'Ka', 15, 'E')
    layer1 = FL(els=[el_p], rho=4.23)
    samp = AnalysisSample(toa=40, volts=[15], layers=[layer1], phimodel='E')
    keV = np.logspace(-2, 2, 100)
    pap_dec = []
    bethe_dec = []
    for V in keV:
        el_p.volt = V
        out = decel_pap(samp, el_p, 'dE/dps')
        pap_dec.append(out[0])
        bethe_dec.append(out[2])
    ax = plt.subplot2grid((1, 1), (0, 0))
    ax.plot(keV, np.asarray(pap_dec), 'b-')
    ax.plot(keV, np.asarray(bethe_dec), 'r-')
    ax.set_xlabel('Accelerating energy [keV]')
    ax.set_ylabel('Average energy loss [keV cm2/g]')
    ax.set_xscale("log")
    #ax.set_yscale("log")
    ax.set_ylim([0, 1e5])
    plt.show()
