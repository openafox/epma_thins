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

from analysis_sample import AnalysisSample as AS
from atomic_element import AtomicElement as AE

def decel_pap(samp, el_p, line=None):
    """Returns Total Trajectory
        or Deceleration Factor if line is provided
        [1] pages 34-6
    """
    # Overvoltage Ratio
    u0 = el_p.e0/el_p.xray

    z = 0.
    mavg = 0.
    javg = 0.
    for lar_i, el_i in samp:
        # p35 eq(6)
        mavg += el_i.z / el_i.mass * el_i.c1
        # p35 eq(7)
        # need to check refs for full explination of where this comes from
        javg += (el_i.c1 * el_i.z / el_i.mass *
            np.log(el_i.z*(10.04 + 8.25*np.exp(-1.0*el_i.z/11.22))/1000.))

    j_mip = np.exp(javg/mavg)  # mean ionization potential of each

    # p35 eq(8)
        p_k = [0.78,
               0.1,
               (-1.*(0.5 - 0.25*j_mip))]

        d_k = [6.6e-6,
               1.12e-5*(1.35 - 0.45*(j_mip**2)),
               2.2e-6/j_mip]
     
    if line is None:
        dr0 = 0.
        for k in range(0,3):
            # p35 eq(9)
            p1 = 1.0 + p_k[k]
            p2 = 1.0 - p_k[k]
            # Total Trajectory [g/cm^2]  [1] page 35
            dr0 += 1/mavg*(j_mip**p2*d_k[k])*(el_p.e0**p1 - el.xray**p1)/p1
            # updated with E0-Eq need to find ref
        return dr0
    elif line[0] == 'K':
        if el_p.z > 30:
            mparam = 0.86
        else:
            mparam = 0.86 + 0.12*np.exp(-1.0*(el_p.z/5.0)**2.0)
    elif line[0] == 'L':
        mparam = 0.82
    elif line[0] == 'M':
        mparam = 0.78
    # if i == 16:  I AM Not sure what this was about   nel was limited to 15??
    #    mparam = 0.86

    v0 = e0/javg
    qe = np.log(u0)/el_p.xray**2/u0**mparam
    df = 0.0
    for k in range(0, 3):
        t_k = 1.0 + p[k] - mparam
        # p36 eq(11) 
        df += [u0/(v0 * mavg)](d_k[k] * (v0/u0)**p[k] *
                   (t_k * u0**t_k * np.log(u0) - u0**t_k + 1.0)/t**2)/qe
        # qe addition need to find ref for this...
    return df

def ppint1(rhoz, L, R):
    """evaluate limits of weighting function integral for PAP model
    10/89, r.a. waldo

    [1] pg54 eq39 with extra factors (/5, /2, /3) where do they come from?
    What is N in the paper???
    """
    ppint1 = rhoz*(rhoz**4/5.0 -
                rhoz**3/2.0*(R + L) +
                rhoz**2/3.0*(R**2 + 4.*R*L + L**2) -
                rhoz*(R**2*L+L**2*R)+(R**2*L**2))
    return ppint1


def ppint2(a_k, b_k, r_k, chi, rhoz):
    """Evaluate limits of x-ray intensity integral.
    [1] p50 eq(35)
    10/89, r.a. waldo
    Arguments:
    a_k -- parabolic branch [1] p37 eq(14)
    b_k -- parabolic branch [1] p37 eq(15)
    r_k -- depth at which the extreem of branch k is situated 
        (Rk=Rinf for first branch, Rk=Rx and Bk=0 for second branch)
    chi -- absorption factor
    rhoz -- mass density
    """
    if chi == 0.0:
        # simplified version to prevent errors
        ppint2 = rhoz*(rhoz**2/3.0 - r_k*rhoz + r_k**2 + b_k/a_k)

    elif chi*rhoz > 88.0:   # Still no Idea where this comes from...
        ppint2 = 0.0
    else:
        ppint2 = np.exp(-chi*rhoz)*((rhoz - r_k)**2 + 2.0*(rhoz - r_k)/chi + 
                                    2.0/chi**2+b_k/a_k)
    return ppint2


def paplimts(ad1,ad2,icase,rc,rx):
    """This subroutine determines the limits of integration for primary and
    secondary x-ray intensities.  Because the PAP model uses two quadratic
    functions to approximate the x-ray depth distribution, the integration
    could be one of three cases;
        1) with limits in phi(a) only
        2)  "     "    "  phi(b)  "
        3)  "   lower limit in phi(a), upper limit in phi(b)

    program completed 11/90 by richard a. waldo

    FIGURE OUT WHERE COMES FROM
    """
    lim = [[0.e-6, rc],[rc, rx]]
    if icase == -1:  # Bulk
        return lim
    elif icase == len(nel) - 1:  # Top Layer
        if rc > ad2:
            lim[0, 1] = ad2
            lim[1, 0] = 1.0
        elif ad2 < rx:
            lim[1, 1] = ad2
    elif icase == 0:  # substrate
        if rc < ad1:
            lim[0, 0] = 1.0
            lim[1, 0] = ad1
        else:
            lim[0, 0] = ad1
    else:  # inner layers
        if rc < ad1:
            lim[0, 0] = 1.0
            lim[1, 0] = ad1
            if ad2 < rx:
                lim[1, 1] = ad2
        else:
            lim[0, 0] = ad1
            if rc > ad2:
                lim[0, 1] = ad2
                lim[1, 0] = 1.0
            else:
                if ad2 < rx:
                    lim[1, 1] = ad2
    if lim[1, 0] > rx:
        lim[1, 0] = 1.0

    return lim


def papint(adelta,icase,da1,da2,db1,drc,drm,drx,c,chiovl,dza,csctheta)
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
    STILL NEED TO UNDERSTAND ND FIX THIS
    """
    db2 = 0.0
    mu = c*csctheta
    rc = drc
    rx = drx
    rm = drm

    # determine the limits of integration of the functions h1, h2
    if icase != -1: # bulk??  need to check goto 20
        if icase == len(nel) - 1: # top layer 
            ad1 = 0.0
            ad2 = adelta[len(nel) - 1]
            #goto 20
        else:
            ad11 = 0.0
            factor = 0.0
            do 10 j=1,icase-1
                ad11=ad11+adelta(j)
                factor = factor + (chiovl[j]*csctheta - mu)*adelta[j]
            if icase == 0:  # substrate
                ad1 = ad11
                ad2 = rx
            else:
                ad1 = ad11
                ad2 = ad11 + adelta[icase]
    # 20
    lim = paplimts(ad1, ad2, icase, rc, rx)

    if lim[0, 0] == 1.0:
        dza1 = 0.0
    else:
        if mu == 0.0:
            dza1 = da1*(ppint2(da1, db1, drm, mu, lim[0, 1]) -
                        ppint2(da1, db1, drm, mu, lim[0, 0]))
        else:
            dza1 = -da1/mu*(ppint2(da1, db1, drm, mu, lim[0, 1]) -
                            ppint2(da1, db1, drm, mu, lim[0, 0]))
    if lim[1, 0] == 1.0:
        dza2 = 0.0
    else:
        if mu == 0.0:
            dza2 = da2*(ppint2(da2, db2, drx, mu, lim[1, 1]) -
                        ppint2(da2, db2, drx, mu, lim[1, 0]))
        else:
          dza2 = -da2/mu*(ppint2(da2, db2, drx, mu, lim[1, 1]) -
                          ppint2(da2, db2, drx, mu, lim[1, 0]))
    za1 = dza1
    za2 = dza2
    dza = dza1 + dza2
    if dza < 0.0:
        dza = 1.e-11
        return dza
    if icase == len(nel) - 1 or icase == -1: # top layer or bulk
        return dza
    dza = dza*np.exp(-1.0*min(factor, 88.0))
    return dza


def papwt(samp, L, R):
    """weighting subroutine for pap parameters
    """
    xnorm = 1.0/(ppint1(R, L, R) - ppint1(0.0, L, R))

    for lar_i, el_i in samp:
        i = sapm.layers.index[lar_i]
        if i == len(samp.layers):  #top layer
            a = 0.0
            b = min([lar_i.thick, R])
        elif i == 0:  #substrate
            a = min([lar_i.depth, R])
            b = R
        else:
            a = min([smap.layers[i + 1].depth, R])
            b = min([lar_i.depth, R])
        wt = xnorm*(ppint1(b, L, R) - ppint1(a, L, R))
        el_i.c1 = el_i.wtfrac * wt
    csum = np.asarray([el_i.c1 for lar_i, el_i in samp])
    el_i.c1 = el_i.c1/csum for lar_i, el_i in samp 

def pap(samp, el_p):
        #i,nel,nels,zs,masses,cnc,e0,xline,line,toa,mac,mode,
        #a,a1,a2,b1,rc,rm,rx,symbol,caller,z)

    # For coding ease
    samp = AS()
    el_p = AE("Mg", "Ka", 15, "s")

    #  initialize variables  
    # cosecant of take off angle
    csctheta = 1.0/np.sin(samp.toa*180/np.pi())
    # Overvoltage Ratio
    u0 = el_p.e0/el_p.xray
    # Path distance or Depth ??
    rt = 7.0e-6*el_p.e0**1.65  # g/cm2  # R in eq 16 in Raw MA '88
    sum1 = 0.0
    # thick0 should have already been run...???
    # get film depths - samp.layer[i].depth
    samp.get_layers_depth()

    # Calculate starting compositions for the weighted-composition-iteration
    # using the same formula in my MAS'88 paper which was used
    # for starting thicknesses.
    if mode == 'F':
        # find element contrabutions
        for lar_i, el_i in samp
            i = samp.layers.index(lar)
            if i = len(samp.layers): # Top Layer 
                el_i.c1 = el_i.wtfrac*(2.0 - lar_i.thick/rt)*lar_i.thick/rt
            elif i == 0:  # substrate  old 7
                el_i.c1 = (el_i.wtfrac*(1.0 - lar_i.depth/rt) *
                         (rt - lar_i.depth)/rt)
                # flmdepth(6)
            else:
                el_i.c1 = (el_i.wtfrac*
                        (((2.0 - lar_i.depth/rt)*lar_i.depth/rt) -
                         ((2.0 - samp.lar_i[i - 1].depth/rt) * 
                          samp.lar_i[i - 1].depth/rt
                          )
                         )
                        )
# Still Checking the above Eqs - Not sure where they are from yet

                # c1[j] = (cnc[j]*(2.-(adelta(layer-1)+adelta(layer))/rt)*
                 #           (adelta(layer)-adelta(layer-1))/rt)
            if el_i.c1 <= 0.0:
                el_i.c1 = 1.e-5
    elif mode == 'B':
        for el_i in samp.layers[0].els:
            el_i.c1 = el_i.wtfrac
    iter1 = 0
    sum1 = np.sum(np.asarray((el_i.c1 for (lar, el_i) in samp)))
    for lar, el_i in samp:
        el_i.c1 = el_i.c1/sum1

    while True:  # 1
        iter += 1
        # Total Trajectory [g/cm^2]  [1] page 35
        dr0 = decel_pap(samp, el_p)

        # ###### Ionization range calc: [1] page 60 appendix 2 ###############
        # Atomic num vars
        z_ln = np.exp(np.sum([np.log(el_i.z)*el_i.c1 for lar_i, el_i in samp]))
        z_mean = np.sum([el_i.z * el_i.c1 for lar_i, el_i in samp])
        # Calc
        q0 = 1.0 - 0.535*np.exp(-(21.0/z_ln)**1.2) - 0.00025*(z_ln/20.0)**3.5
        beta = 40.0/z_mean
        operand = min([(u0-1.0)/beta,88.0])  # why limit to 88 ? fromGRMFilm
        # Emperical factor adjusting the Range of Ionization
        Q_c = q0 + (1.0 - q0)*np.exp(-operand)
        # factor
        h_c = z_mean**0.45  # adelta
        # D factor accounting for dispersion af traj at low u0
        d_lu0 = 1.0 + 1.0/(u0**h_c)
        # Range of Ionization
        drx = Q_c*d_lu0*dr0
        rx = drx
        # ####################################################################

        if mode != 'B':
            # if film do the pap weighting functs  NEED TO LOOK AT!!!!
            if np.abs(rx*1.e6 - rt*1.e6) > 0.02:
                c1 = papwt(h_c, cnc, nel, nels, -0.4*rx, rx)
                #starting weighting
                # pg54 is it L = -0.4R(paper) or L= -0.4Rx(here)?
                rt = rx
                # goto 1
            else:
                break
            # get bulk concentrations
            c1 = papwt(samp, -0.4*rx, 0.5*rx)
            #near surface...???
            # pg54 is it L = -0.4R(paper) or L= -0.4Rx(here)?
        else:
            break

    # [1] page 59 appendix 1 (Backscatter Loss Factor) ##############
    # z mean averaging rule
    zp = np.sum([el_i.z**0.5*el_i.c1 for lar_i, el_i in samp])**2
    # Coefficient
    eta = 0.00175*zp + 0.37*(1.0 - np.exp(-0.015*(zp**1.3)))
    # BS coefficient (Er/E0)
    wavg = 0.595 + eta/3.7 + eta**4.55

    q_gu = (2.0*wavg - 1.0)/(1.0 - wavg)
    
    ju = 1.0 + u0*(np.log(u0) - 1.0)
    
    gu = ((u0 - 1.0 - ((1.0 - (1.0/u0**(q_gu + 1.0)))/(1.0 + q_gu)))/
          ((2.0 + q_gu)*ju))
    # Backscatter Loss Factor
    rb = 1.0 - eta*wavg*(1.0 - gu)
    # ###############################################################

    # Surface Ionization ([1] page 60 appendix 2)
    dphi0 = 1.0 + 3.3*(1.0 - (1.0 /u0**(2.0 - 2.3*eta)))*eta**1.2

    # get init concentrations
    if mode == 'F':
        c1 = papwt(samp, -0.6*rx, 0.7*rx)
        # Leads to mean atomic number in 1/S (decel) calc??
        # pg54 is it L = -0.6R(paper) or L= -0.6Rx(here)?

    # Ionization Crossection
    df = rb*decel_pap(samp, el_p, el_p.line, rb)

    if mode == 'F':
        z_mean = (z_mean + zp)/2.0   # Not sure why we do this...

    # Depth of Maximum Phi(rho z) ([1] page 61 appendix 2) ###################
    g1z = 0.11 + 0.41*np.exp(-(z_mean/12.75)**0.75)

    g2z = 1.0
    t1 = -(u0-1.0)**0.35/1.19
    if t1 > -88.:
        g2z = 1.0 - np.exp(t1)

    g3z = 1.0
    t2 = -(u0 - 0.5)*z_mean**0.4/4.0
    if t2 > -88.0:
        g3z = 1.0 - np.exp(t2)

    drm = g1z*g2z*g3z*drx #drx = range of ionization (above)
    # ######################################################################

    # p38 eq(17)
    delt = ((drx - drm) * (df - dphi0 * drx/3.0) *
            ((drx - drm) * df - dphi0 * drx * (drm + drx/3)))
    if delt < 0.0:
        # Modifications for small overvoltages (p61 Appendix 3)
        rmt = drm
        drm = drx * (df - dphi0 * drx/3.0)/(df + dphi0 * drx)
        rm = drm
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
                   % (el, line, rmt*1.0e6, rm*1.0e6, callerb))

    # Root (rho z where two parabolas meet) p38 eq(17)
    drc = 1.5*((df - dphi0 * drx/3.0)/dphi0 - 
               np.sqrt(delt)/(dphi0 * (drx - drm)))
    # Coefficients of the parobolic branches (p38 eq(18)
    da1 = dphi0/(drm*(drc - drx*(drc/drm - 1.0)))
    da2 = da1*(drc - drm)/(drc - drx)
    db1 = dphi0 - da1*drm*drm

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

    papint(adelta,layer,da1,da2,db1,drc,drm,drx,
     &chi,chiovl,dza1,csctheta)
    za(2) = df
    a1 = da1
    a2 = da2
    b1 = db1
    rm = drm
    rx = drx
    rc = drc
    za(1) = dza1
         
         
if __name__ == '__main__':
    pass

