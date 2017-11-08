#!/usr/bin/env python
"""-----------------------------CHARFLR-------------------------------

This program controls the calculation of the CHARacteristic
x-ray FLuoRescence correction for thin film and bulk systems.
program completed on 6/88 by richard a. waldo

Keyword arguments:
A -- applena,nel,nels,mode,zx,ax,symb,e0,cnc,delta,fchar,phipar,toa,mac,exciter,fmac)
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import numpy as np
from layr import layr

def charflr(na, nel, nels, mode, zx, ax, symb, e0, cnc, delta, fchar, phipar,
        toa, mac, exciter, fmac):
      # real zx(15),ax(15),cnc(15),mac(15,15),delta(7),
     # &fcp(7),fcs(7),delta2(9),za(2),exciter(200,3),fmac(15,200)
      # integer nel(7)
      # jlogical ff
# ff overrides R(b/2a) warning message in subroutine prozbeta
      # character*1 phipar,mode
      # character*7 symb(15),symbol
      # character*3 lines(12)
    fcs = []
    lines = ['Kb ', 'Ka ', 'Lg2', 'Lb3', 'Lb4', 'Lg1',
             'Lb1', 'Lb2', 'La1', 'Mg ', 'Mb ', 'Ma ']
    fchar = 0.0
    nnsum = 0
    csctheta = 1/np.sin(toa/57.29578)
    ff = True 
    delta[6] = 1.0
    rx = 1.5*6.5*(e0**1.7)/1.0e6
    temp = 0.0
    delta2[0] = 0.0
    for kk in range(1, 6):
        temp = temp + delta[kk - 1]
        delta2[kk] = temp
    layera = layr(nel, na)
    for ii in range(0, 6):
        fcs.append(0.0)
    for ik in range(1, nels):
        layerc = layr(nel, ik)
        fcs[layerc] = fcs[layerc] + mac[na, ik]*cnc[ik]*csctheta
    nn1 = int(exciter[184 + na,1])
    nn2 = int(exciter[185 + na,1])
    if nn2 - nn1 = 0: 
        return
    for nlb in range(nn1+1, nn2):
        j = int((exciter(nlb, 3)-1.0)/12) + 1
        lb = int(exciter(nlb, 3)-12.0*(j - 1))
        xline = exciter(nlb, 2)
        layerb = layr(nel, j)
        leb = nedge[lb]
        symbol[1:2] = symb(j)(1:2)
        symbol[3:5] = lines(lb)
        ###  hereh  ###
        if phipar == 'E':
            call pap(delta,j,nel,nels,zx,ax,cnc,e0,xline,lb,
     &    toa,mac,mode,za,a1,a2,b1,rc2,rm,rx,symbol,mode,z)
        else
          call phirzeq(nel,nels,j,delta,alpha,beta,gamma,phi0
     &    ,zx,ax,cnc,e0,xline,phipar,lb,ff,z)
        endif
        do 410 ii=1,7
410       fcp(ii)=0.
        do 420 ik=1,nels
          layerc=layr(nel,ik)
420       fcp(layerc)=fcp(layerc)+cnc(ik)*fmac(ik,nlb)
        if (mode.eq.'B') then
          delta2(8)=0.
          delta2(9)=rx
          else if (mode.eq.'F') then
            if (layerb.eq.7) then
              delta2(8)=delta2(7)
              delta2(9)=rx
            else
              delta2(8)=delta2(layerb)
              delta2(9)=delta2(layerb+1)
            endif
            delta2(9)=amin1(delta2(9),rx)
            if (delta2(8).ge.delta2(9)) goto 5000
          endif
        if (phipar.eq.'E') then
          call papfluor(delta2,fcp,fcs(layera),layera,layerb,
     &    ffact,a1,a2,b1,rc2,rm,rx,mode)
        else
          call tripint(delta2,fcp,fcs(layera),layera,layerb,
     &    ffact,alpha,beta,gamma,phi0,mode)
        endif
        ffact1=ffact
        factor=0.
        if (mode.eq.'B') goto 771
        do 770 kk=2,7
          factor=factor+fcs(kk-1)*delta(kk-1)
          if (layera.eq.kk) ffact=ffact*exp(-1.*amin1(factor,88.))
770     continue
771     continue
        ffact=0.5*cnc(na)*cnc(j)*ffact*exciter(nlb,1)
        fchar=fchar+ffact
5000  continue
      return
10000 end
if __name__ == '__main__':
    pass

