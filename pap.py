#!/usr/bin/env python
"""
Calculates x-ray intensities using PAP model.

[1]J.-L. Pouchou and F. Pichoir, “Quantitative analysis of homogeneous or
stratified microvolumes applying the model ‘PAP’”
in Electron probe quantitation, Springer, 1991, pp. 31–75.

Vars of note:
c1 --  ??  weight concentration of analyzed element
samp.layers[i].thick -- adelta --
cnc -- weightpercent
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import numpy as np
from layrelem import layrelem

from atomic_element import AtomicElement as AtEl
from film_layer import FilmLayer as FL
from analysis_sample import AnalysisSample as AS

# for coding purposes
ti = AtEl('Ti', 'Ka', 15, 's')
layer1 = FL(els=[Si, o], rho=2.65)
layer2 = FL(els=[ti, o], rho=4.23)
samp = AS(toa=40, volts=[15], layers=[layer1, layer2], phimodel='E')
####################


def decel_pap(samp, el_p, line=None):
    """Returns Total Trajectory
        or Deceleration Factor if line is provided
        el_p -- primary element ie element of intrest of excited element
        [1] pages 34-6
    """
    # Overvoltage Ratio
    u0 = el_p.e0/el_p.xray

    z = 0.
    mavg = 0.
    javg = 0.
    for lar_i, el_i in samp:
        # [1]p35e6
        mavg += el_i.z / el_i.mass * el_i.c1
        # [1]p35e7
        # need to check refs for full explination of where this comes from
        javg += (el_i.c1 * el_i.z / el_i.mass *
                 np.log(el_i.z*(10.04 + 8.25*np.exp(-1.0*el_i.z/11.22))/1000.))

    j_mip = np.exp(javg/mavg)  # mean ionization potential of each [1]p35e6

    # [1]p35e8
    p_k = [0.78,
           0.1,
           (-1.*(0.5 - 0.25*j_mip))]

    d_k = [6.6e-6,
           1.12e-5*(1.35 - 0.45*(j_mip**2)),
           2.2e-6/j_mip]

    if line is None:
        dr0 = 0.
        for k in range(0, 3):
            # [1]p35e9
            p1 = 1.0 + p_k[k]
            p2 = 1.0 - p_k[k]
            # Total Trajectory [g/cm^2]  [1]p35e9
            dr0 += 1/mavg*(j_mip**p2*d_k[k])*(el_p.e0**p1 - el_p.xray**p1)/p1
            # updated with E0-Eq need to find ref
        return dr0
    # Ionization crosssection [1]p36e10-11
    # mparam from Hutchins (get Ref)
    elif line[0] == 'K':
        mparam = 0.86 + 0.12*np.exp(-1.0*(el_p.z/5.0)**2.0)
    elif line[0] == 'L':
        mparam = 0.82
    elif line[0] == 'M':
        mparam = 0.78

    v0 = el_p.e0/javg
    qe = np.log(u0)/(el_p.xray**2*u0**mparam)  # [1]p36e10
    df = 0.0
    for k in range(0, 3):
        t_k = 1.0 + p[k] - mparam
        # [1]p36e11   (1/S)*qe [1]p36e13
        df += ([u0/(v0 * mavg)](d_k[k] * (v0/u0)**p[k] *
               (t_k * u0**t_k * np.log(u0) - u0**t_k + 1.0)/t**2)/qe)
    return df


def ppint1(rhoz, L, R):
    """evaluate limits of weighting function integral for PAP model
    10/89, r.a. waldo

    [1]p54e39 with extra factors (/5, /2, /3) where do they come from?
    What is N in the paper???
    """
    ppint1 = rhoz*(rhoz**4/5.0 -
                   rhoz**3/2.0*(R + L) +
                   rhoz**2/3.0*(R**2 + 4.*R*L + L**2) -
                   rhoz*(R**2*L+L**2*R)+(R**2*L**2))
    return ppint1


def ppint2(a_k, b_k, r_k, chi, rhoz):
    """Evaluate limits of x-ray intensity integral.
    [1]p50e35
    10/89, r.a. waldo
    Arguments:
    a_k -- parabolic branch [1]p37e14
    b_k -- parabolic branch [1]p37e15
    r_k -- depth at which the extreem of branch k is situated
        (Rk=Rinf for first branch, Rk=Rx and Bk=0 for second branch)
    chi -- absorption factor
    rhoz -- mass density
    """
    if chi == 0.0:
        # simplified version to prevent errors
        ppint2 = a_k*rhoz*(rhoz**2/3.0 - r_k*rhoz + r_k**2 + b_k/a_k)

    elif chi*rhoz > 88.0:   # rounding to 0 exp->0  why 88?? we may never know
        ppint2 = 0.0
    else:
        ppint2 = (-(a_k/chi)*np.exp(-chi*rhoz)*((rhoz - r_k)**2 +
                  2.0*(rhoz - r_k)/chi + 2.0/chi**2+b_k/a_k))
    return ppint2


def papint(samp, lar_p, a_1, a_2, b_1, r_c, r_m, r_x, c, chiovl, dza,
           csctheta):
    """This program integrates the phi(rz) curve for PAP model.  4 cases
    need to be handled.  case 1) the phi(rz) curve of a bulk
    specimen (for example when the pure element standard x-ray
    intensity is needed to calculate a k-ratio.  case 2) the
    curve for a thin film of thickness delta.  case 3)  the curve
    for a substrate under a film of thickness delta, the curve is
    integrated from delta to infinity.  case 4) the curve for a buried
    film of thickness delta1.  for thin films, weighted
    average values of the phi(rz) parameters are used (see
    subprogram to see how the parameter vlues are weighted).

    program completed 10/89 by richard a. waldo
    rc -- intersection of 2 parabolic curves
    r_x -- max depth (rhoz) aka ionization range
    r_m -- position (rhoz) at max height (phi)
    """

    db2 = 0.0
    mu = lar_p.chi * csctheta

    # MAY NEED IF BULK CASE!!!
    # [[0.e-6, r_c],[r_c, r_x]]

    # upper integration bound (bottom of layer, top of next bellow)
    u_bound = samp.layers[l_index - 1].depth
    # lower integration bound (top of layer)
    l_bound = lar_p.depth

    if u_bound > r_x:
        # don't integrate beyond the ionaziation range
        u_bound = r_x

    if u_bound < r_c:
        # evaluate only first int
        dza1 = (ppint2(a_1, b_1, r_m, mu, u_bound) -
                ppint2(a_1, b_1, r_m, mu, l_bound))
        dza2 = 0.0
    elif u_bound > r_c and l_bound < r_c:
        # evaluate both
        # first with l_bound to r_c
        dza1 = (ppint2(a_1, b_1, r_m, mu, r_c) -
                ppint2(a_1, b_1, r_m, mu, l_bound))
        # second with r_c to u_bound
        dza2 = (ppint2(a_2, db2, r_m, mu, u_bound) -
                ppint2(a_2, db2, r_m, mu, r_c))

    elif l_bound > r_c:
        # only eval upper int
        dza1 = 0.0
        dza2 = (ppint2(a_2, db2, r_m, mu, u_bound) -
                ppint2(a_2, db2, r_m, mu, l_bound))
    else:
        pass
        # Alert Error

    za1 = dza1
    za2 = dza2
    dza = dza1 + dza2
    if dza < 0.0:
        dza = 1.e-11
        return dza

    # Get T_A factor for [1]p50e35 (its on p 49)
    factor = 0.0
    for j in range(len(samp.layers), samp.layers.index(lar_p)):
        factor += (samp.layers[j].thick*(samp.layers[j].chiovl * csctheta -
                   lar_p.chi * csctheta))

    if icase == len(nel) - 1 or icase == -1:  # top layer or bulk
        return dza
    dza = dza*np.exp(-1.0*min(factor, 88.0))
    return dza


def papwt(samp1, L, R):  # remove the 1 after coding
    """weighting subroutine for pap parameters
    """
    xnorm = 1.0/(ppint1(R, L, R) - ppint1(0.0, L, R))

    for lar_i, el_i in samp:
        i = samp.layers.index[lar_i]
        if i == len(samp.layers):  # top layer
            a = 0.0
            b = min([lar_i.thick, R])
        elif i == 0:  # substrate
            a = min([lar_i.depth, R])
            b = R
        else:
            a = min([samp.layers[i + 1].depth, R])
            b = min([lar_i.depth, R])
        wt = xnorm*(ppint1(b, L, R) - ppint1(a, L, R))
        el_i.c1 = el_i.wtfrac * wt

    csum = np.asarray([el_i.c1 for lar_i, el_i in samp])
    el_i.c1 = el_i.c1/csum for lar_i, el_i in samp


def pap(samp1, el_p):
    # initialize variables

    # cosecant of take off angle
    csctheta = 1.0/np.sin(samp.toa*180/np.pi())
    u0 = el_p.e0/el_p.xray
    # Ionization Depth (Max depth)
    rt = 7.0e-6*el_p.e0**1.65  # [g/cm2] # R in eq 16 [RAW-MAS-88]
    # would it be more accurate to use the Eo-Ec eq7 in [RAW-MAS-88]
    sum1 = 0.0

    # thick0 should have already been run...???

    # get film depths - samp.layer[i].depth
    samp.get_layers_depth()

    # Calculate starting compositions [RAW-MAS-88]
    if mode == 'F':
        # find element contrabutions
        for lar_i, el_i in samp:
            i = samp.layers.index(lar)
            if i = len(samp.layers):  # Top Layer
                # DOCS [RAW - starts] eq2
                el_i.ci = el_i.ci*(2.0*lar_i.thick/rt - lar_i.thick**2/rt**2)
            elif i == 0:  # substrate
                # DOCS [RAW - starts] eq4 R subed for y
                el_i.ci = el_i.ci*(1.0 - 2*lar_i.depth/rt +
                                   (lar_i.depth/rt)**2)
            else:  # Burried Layer
                # DOCS [RAW - starts] eq4
                el_i.ci = (el_i.ci *
                           (((2.0 - samp.layers[i - 1].depth/rt) *
                             samp.layers[i - 1].depth/rt) -
                            ((2.0 - samp.lar_i.depth/rt) *
                             samp.lar_i.depth/rt)
                            )
                           )
            if el_i.c1 <= 0.0:
                el_i.c1 = 1.e-5
            # Log ci
            el_i.clog.append[el_i.ci, 'start comp']
    # elif bulk no change c=wtfrac <- both stored at same place...?

    iter1 = 0
    # Normalize C1 eq13 (eqs above are c1=c1*k??)
    sum1 = np.sum(np.asarray((el_i.c1 for (lar, el_i) in samp)))
    for lar, el_i in samp:
        el_i.c1 = el_i.c1/sum1

    while True:  # 1
        iter += 1
        # Total Trajectory [g/cm^2]  [1]p35
        dr0 = decel_pap(samp, el_p)

        # ###### Ionization range calc: [1]p60 appendix 2 ###############
        # Atomic num vars
        z_ln = np.exp(np.sum([np.log(el_i.z)*el_i.c1 for lar_i, el_i in samp]))
        z_mean = np.sum([el_i.z * el_i.c1 for lar_i, el_i in samp])
        # Calc
        q0 = 1.0 - 0.535*np.exp(-(21.0/z_ln)**1.2) - 0.00025*(z_ln/20.0)**3.5
        beta = 40.0/z_mean
        operand = min([(u0-1.0)/beta, 88.0])  # why limit to 88 ? fromGRMFilm
        # Emperical factor adjusting the Range of Ionization
        Q_c = q0 + (1.0 - q0)*np.exp(-operand)
        # factor
        h_c = z_mean**0.45
        # D factor accounting for dispersion af traj at low u0
        d_lu0 = 1.0 + 1.0/(u0**h_c)
        # Range of Ionization
        r_x = Q_c*d_lu0*dr0
        # ####################################################################

        if mode != 'B':
            # in GMRFilm weighing the initial C (C_m-1, cnc) not new Ci
            # and not the new one (C_m, c1)???????????
            # Changed here to el_ci because it make more sense to use it but
            #maybe I'm wrong?????

            if np.abs(r_x*1.e6 - rt*1.e6) > 0.02:
                # stop if diffrence is less than 1%??  is 0.02 = 1%?
                c1 = papwt(samp, -0.4*r_x, r_x)
                # starting weighting
                # [1]p54
                rt = r_x
            else:
                c1 = papwt(samp, -0.4*r_x, 0.5*r_x)
                #near surface...???
                # pg54 is it L = -0.4R(paper) or L= -0.4Rx(here)?
                # above it did not matter bc R=Rx but here R=Rx/2
                break  # stop loop
        else:
            break

    # [1]p59 appendix 1 (Backscatter Loss Factor) ##############
    # z mean averaging rule
    zb = np.sum([el_i.z**0.5*el_i.c1 for lar_i, el_i in samp])**2
    # Coefficient
    eta = 0.00175*zb + 0.37*(1.0 - np.exp(-0.015*(zb**1.3)))
    # BS coefficient (Er/E0)
    wavg = 0.595 + eta/3.7 + eta**4.55

    q_gu = (2.0*wavg - 1.0)/(1.0 - wavg)

    ju = 1.0 + u0*(np.log(u0) - 1.0)

    gu = ((u0 - 1.0 - ((1.0 - (1.0/u0**(q_gu + 1.0)))/(1.0 + q_gu))) /
          ((2.0 + q_gu)*ju))
    # Backscatter Loss Factor
    rb = 1.0 - eta*wavg*(1.0 - gu)
    # ###############################################################

    # Surface Ionization ([1]p60 appendix 2)
    dphi0 = 1.0 + 3.3*(1.0 - (1.0/u0**(2.0 - 2.3*eta)))*eta**1.2

    if mode == 'F':
        c1 = papwt(samp, -0.6*r_x, 0.7*r_x)
        # Leads to mean atomic number in 1/S (decel) calc??
        # pg54 is it L = -0.6R(paper) or L= -0.6Rx(here)?

    # Ionization Crossection
    df = rb*decel_pap(samp, el_p, el_p.line)

    if mode == 'F':
        z_mean = (z_mean + zb)/2.0   # Not sure why RAW did this...
        # I currently have not found it in any of the PAPers

    # Depth of Maximum Phi(rho z) ([1] page 61 appendix 2) ###################
    g1z = 0.11 + 0.41*np.exp(-(z_mean/12.75)**0.75)
    g2z = 1.0 - np.exp(-(u0-1.0)**0.35/1.19)
    g3z = 1.0 - np.exp(-(u0 - 0.5)*z_mean**0.4/4.0)

    r_m = g1z*g2z*g3z*r_x
    # ######################################################################

    # [1]p38e17
    delt = ((r_x - r_m) * (df - dphi0 * r_x/3.0) *
            ((r_x - r_m) * df - dphi0 * r_x * (r_m + r_x/3)))
    if delt < 0.0:
        # Modifications for small overvoltages ([1]p61 Appendix 3)
        rmt = r_m
        r_m = r_x * (df - dphi0 * r_x/3.0)/(df + dphi0 * r_x)
        delt = 0.0
        callb = {'E': 'pure element standard',
                 'F': 'fluorescence, film',
                 'C': 'compound standard',
                 'B': 'fluorescence, bulk (std.?)',
                 'I': 'main iteration procedure'}
        callerb = callb(caller)
        #  Warning statements suppressed when calling
        #  routine is either bulk or film fluorescence

        if caller == 'E' or caller == 'C' or caller == 'I':
            print ('!!!Warning!!! overvoltage ratio is too low for %s %s\n'
                   'Rm lowered from %7.2f to %7.2f to calculate parameters;'
                   '\nCalling routine:%s'
                   % (el, line, rmt*1.0e6, r_m*1.0e6, callerb))

    # Root (rho z where two parabolas meet) [1]p38e17
    r_c = 1.5*((df - dphi0 * r_x/3.0)/dphi0 -
               np.sqrt(delt)/(dphi0 * (r_x - r_m)))
    # Coefficients of the parobolic branches p38e18
    a_1 = dphi0/(r_m*(r_c - r_x*(r_c/r_m - 1.0)))
    a_2 = a_1*(r_c - r_m)/(r_c - r_x)
    b_1 = dphi0 - a_1*r_m*r_m

    """  I have no idea what this was for....
    if i != 16: # goto 451
        ### Need to figure out this i1, i2 mess
        i1 = 1
        i2 = nel[len(nel) - 1]
        layer = layrelem(nel,i)
        if caller == 'F' or caller == 'B':
            chi = 0.0
            chiovl = np.zeros(len(nel))
            # do 447 n=1,7
            # chiovl(n)=0.
        elif mode != 'B':  # goto 449
            if layer >= 2:
                chiolv = chiov(mac, i, nel, layer, cnc)
            for kk in range(0,
                if (layer != kk):
                    i1 = i1 + nel[kk]
                    i2 = i2 + nel[kk+1]
        # 449
        chi = 0.0
        for j in range(i1, i2):
            chi = chi + cnc[j]*mac(i,j)
        # 451
        if mode = 'B':
            layer = 0 # 8  #???  in orriginal below sub???
        """

    papint(adelta, layer, a_1, a_2, b_1, r_c, r_m, r_x,
           chi, chiovl, dza1, csctheta)
    za(2) = df
    za(1) = dza1

    return (a_1, b_1, a_2, r_m, r_c, r_x)


if __name__ == '__main__':
    from atomic_element import AtomicElement as AtEl
    from film_layer import FilmLayer as FL
    Si = AtEl('Si', 'Ka', 15, 's')
    o = AtEl('O', 'Ka', 15, 's')
    ti = AtEl('Ti', 'Ka', 15, 's')
    layer1 = FL(els=[Si, o], rho=2.65)
    layer2 = FL(els=[ti, o], rho=4.23)
    samp = AnalysisSample(toa=40, volts=[15], layers=[layer1, layer2],
                          phimodel='E')
    print 'defs done'
    for lay, el in samp:
        print ('Element:', el.name, ' Macs', el.mac, 'layer',
               samp.layers.index(lay))
    print ""
    for lay, el in samp:
        print 'Element:', el.name, 'Density', lay.rho
    # Need to add test of calc_thick0...
