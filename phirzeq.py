#!/usr/bin/env python
"""-------------------------------PHIRZEQ---------------------------------

calculates the four phi(rz) parameters alpha, beta, gamma
and phi(0) for each layer in a thin film system.  several
models are included (see main program for list).  the variable
phipar contains the i.d. of the particular model used.
program completed 4/87  by richard a. waldo

updated 10/88 for Bastin's nbs workshop model (1988), r.a. waldo
NOTE: see also Scanning vol.12, 1990. p.225

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import numpy as np
from scipy.special import erfc
from scale import scale


def phirzeq(samp, el_p, lay_p, ff):
    #what is ff - boolan
    e0 = el_p.volt
    xline = el_p.xray/1000
    u0 = e0/xline
    alpha = 1.e6/(3.*e0**1.7)
    javg = 0.
    bsf = 0.
    zp = 0.
    for k in range(0, 5):
        ajp = 0.
        z = 0.
        a = 0.
        za = 0.
        talpha = alpha
        scale(samp, talpha)
        for lay_j, el_j in samp:
            if el_j.c1 == 0.:
                    el_j.c1 = 1.e-5
            if samp.phimodel == 'B':
                "This is the Refrence for this"
                el_j.aj = (9.29*el_j.z*(1.+1.287/el_j.z**0.667))/1000.
            elif samp.phimodel == 'C':
                "This is the Refrence for this"
                el_j.aj = el_j.z*(10.04+8.25*np.exp(-1.*el_j.z/11.22))/1000.
            elif samp.phimodel == 'P':
                "This is the Refrence for this"
                ajp = ajp + lay_j.lwt*el_j.c1*.0115*el_j.z

            za = za+(el_j.z/el_j.mass)*lay_j.lwt*el_j.c1
            z = z+el_j.z*lay_j.lwt*el_j.c1
            a = a+el_j.mass*lay_j.lwt*el_j.c1
        alpha = 0.
        for lay_i, el_i in samp:
            if samp.phimodel == 'B':
                "This is the Refrence for this"
                alpha += ((((175000./(e0**1.25*(u0-1.)**0.55)) *
                            ((np.log(1.166*e0/el_i.aj)/xline)**0.5))**2) *
                          (lay_i.lwt*el_i.c1*el_i.z/el_i.mass))
            elif samp.phimodel == 'C':
                alpha += (1./((216140.*el_i.z**1.163 /
                               (e0**1.25*el_i.mass*(u0-1.)**0.5)) *
                              ((np.log(1.166*e0/el_i.aj)/xline)**0.5)) *
                          lay_i.lwt*el_i.c1*el_i.z/el_i.mass)
        if samp.phimodel == 'B':
            "This is the Refrence for this"
            alpha = np.sqrt(a*alpha/z)
        elif samp.phimodel == 'C':
            "This is the Refrence for this"
            alpha = za/alpha
        elif samp.phimodel == 'P':
            "This is the Refrence for this"
            alpha = ((415000./e0**0.75)*((z-1.3)/z)*(za**0.5) *
                     (np.sqrt(za*np.log((1.166/ajp*(e0+xline)/2.)) /
                              (e0**2-xline**2))))
        if samp.layers[-1].thick == 0:  # if top layer thickness = 0
            break
        if (abs(alpha-talpha) < 1.):
            break
    if samp.phimodel == 'B':
        """calculate according to bastin scanning'86 paper"""
        for lay_i, el_i in samp:
            eta = (-52.3791 + 150.48371*el_i.z -
                   1.67373*el_i.z**2 + 0.00716*el_i.z**3)/10000.
            geta = (-1112.8 + 30.289*el_i.z - 0.15498*el_i.z**2)/10000.
            eta = eta*(1.0 + geta*np.log(e0/20.))
            bsf = bsf + lay_i.lwt2*el_i.c1*eta
        gu = -0.59299 + 21.55329/u0 - 30.55248/u0**2 + 9.59218/u0**3
        iu = 3.43378 - 10.7872/u0 + 10.97628/u0**2 - 3.62286/u0**3
        phi0 = 1. + (bsf/(1. + bsf))*(iu + gu*np.log(1. + bsf))
#
        if (u0 <= 3.):
            gamma = 1. + (u0 - 1)/(0.3384 + 0.4742*(u0 - 1))
        if (u0 > 3.):
            gamma = (5.*3.14159*(u0 + 1.)/u0/(np.log(u0 + 1.0)) *
                     (np.log(u0 + 1.) - 5. + 5.*(u0 + 1.)**(-0.2)))
        xn = z/(0.4765 + 0.5473*z)
        beta = alpha*z**xn/a
        return

    if samp.phimodel == 'P':
        """calculate according to packwood"""
        for lay_i, el_i in samp:
            eta = (-52.3791 + 150.48371*el_i.z +
                   (-1.67373)*el_i.z**2 + 0.00716*el_i.z**3)/10000.
            bsf += lay_i.lwt2*el_i.c1*eta
        phi0 = 1. + 2.7*bsf*(1. - np.exp((1-u0)/2.))
        gamma = 31.4159*u0/(u0 - 1.)*(1. + (10./(np.log(u0))*(u0**(-0.1)-1.)))
        beta = 0.4*alpha*z**0.6

    if samp.phimodel == 'C':
        """this section calculates the bastin scanning'90 model parameter
        values."""
        #if last element mparam = 0.86??
        if (el_p.line[0] == 'L'):
            mparam = 0.82
        if (el_p.line[0] == 'M'):
            mparam = 0.78
        if (el_p.line[0] == 'K'):
            if (el_p.z > 30):
                mparam = 0.86
            else:
                mparam = 0.86 + 0.12*np.exp(-1.*(el_p.z/5.)**2.)
        for lay_i, el_i in samp:
            zp += el_i.z**0.5*lay_i.lwt2*el_i.c1
            javg += lay_i.lwt*el_i.c1*el_i.z/el_i.mass*np.log(el_i.mass)
        zp = zp**2
        javg = javg/za
        javg = np.exp(javg)
        pk = [0.78, 0.1, -1.*(0.5-0.25*javg)]
        dk = [6.6e-6, 1.12e-5*(1.35 - 0.45*(javg**2.)), 2.2e-6/javg]
        # calculate phi(0), surface ionization; rb, backscatter factor;
        # and f, the integral of the intensity
        eta = 0.00175*zp + 0.37*(1. - np.exp(-1.*(0.015*(zp**1.3))))
        wavg = 0.595 + eta/3.7 + eta**4.55
        alph = (2.*wavg-1.)/(1. - wavg)
        ju = 1. + u0*(np.log(u0) - 1.)
        gu = ((u0 - 1.0 - ((1.0 - (1.0/u0**(alph+1.0)))/(1.0 + alph))) /
              ((2.0 + alph)*ju))
        gamm = 2.0 - 2.3*eta
        # rb is backscatter factor
        rb = 1.0 - eta*wavg*(1. - gu)
        phi0 = 1. + 3.3*(1.0 - (1./u0**gamm))*eta**1.2
        # qe is ionization cross section at eO
        qe = np.log(u0)/xline**2/u0**mparam
        v0 = e0/javg
        t7 = rb*(u0/v0)/za
        fk = 0
        for k in range(0, 3):
            t = 1.0 + pk[k]-mparam
            t1 = u0**(t)
            fk += t7*(dk[k]*(v0/u0)**pk[k]*(t*t1*np.log(u0)-t1 + 1.)/t**2)/qe
        if (u0 <= 6.0):
            gamma = (3.98352*u0**(-0.0516861)*(1.276233-u0 **
                     (-1.25558*z**(-0.1424549))))
        elif (u0 > 6.0):
            gamma = 2.814333*u0**(0.262702*z**(-0.1614454))
        if (xline < 0.7):
            gamma = gamma*xline/(-.0418780 + 1.05975*xline)
            (beta, flag) = bastbeta(alpha, gamma, phi0, fk)
        if flag:
            if not el_p.flag:
                el_p.flag = True
                if (not ff):
                    print(' !!! 0>R(b/2a)>1 for element Z = %d, xline = %s !!!'
                          '\n!!! Overvoltage ratio may be too low !!!'
                          % (el_p.z, el_p.line))
    return alpha, beta, gamma, phi0


def bastbeta(a, g, p, f):
    """this program finds the proper beta which fits the
    conditions required for bastin's nbs workshop (1988) phi(rz) model
    program completed 10/88 r.a. waldo

    NOTE:see Scanning, vol. 12, 1990, p.225 for this model.

    Note: bastbeta.tmp : see line change below; 5-21-90 r.a.w."""
    flag = False
    x = (g - 2.0*a*f/np.sqrt(3.1415926))/(g - p)
    if x > 1.0 or x < 0.0:
        # Note: bastbeta.ftn (Bastin's formulation) had the following two lines
        # g=(g+p)/2.
        # p=g
        # could this change be more physically realistic
        g = p
        a = (g + p)/np.sqrt(3.1415926)/4.0/f
        x = 0.5
        flag = True
    y1 = rb2a(x)
    beta = y1*(2*a)
    while True:
        beta1 = beta
        y2 = erfc(y1)/x*y1
        beta = y2*(2*a)
        if abs(beta1 - beta) < 1.0:
            return (beta, flag)
        y1 = y2


def rb2a(x):
    """bastins's nbs workshop (1988) phi(rz) model
    NOTE: see also Scanning, Vol.12, 1990, p225."""
    if x < 0.03165:
        b = (1.0 - 4.894396*x)/(1.341313*x)
    elif x < 0.056:
        b = (1.0 - 2.749786*x)/(1.447465*x)
    elif x < 0.102:
        b = (1.0 - 1.043744*x)/(1.604820*x)
    elif x < 0.306:
        b = (1.0 - 0.5379956*x)/(1.685638*x)
    elif x < 0.57:
        b = 4.852357*np.exp(-3.680818*x)
    elif x < 0.7:
        b = 5.909606*np.exp(-4.015891*x)
    elif x < 0.8:
        b = 13.4381*np.exp(-5.180503*x)
    elif x < 0.9:
        b = 1.122405 - 1.141942*x
    else:
        b = 0.9628832 - .9642440*x
    return b

if __name__ == '__main__':
    from atomic_element import AtomicElement as AtEl
    from film_layer import FilmLayer as FL
    from analysis_sample import AnalysisSample
    el_p = AtEl('Si', 'Ka', 15, 'E')
    el2 = AtEl('O', 'Ka', 15, 'S')
    el3 = AtEl('Ti', 'Ka', 15, 'E')
    el4 = AtEl('O', 'Ka', 15, 'E')
    layer1 = FL(els=[el_p, el2], rho=2.65)
    layer2 = FL(els=[el3, el4], rho=4.23)
    samp = AnalysisSample(toa=40, volts=[15], layers=[layer1, layer2],
                          phimodel='E')
    phirzeq(samp, el_p, layer1, False)
