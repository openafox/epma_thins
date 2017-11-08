#!/usr/bin/env python
"""-----------------------------CHARFDTA-------------------------------
CHARacteristic Fluorescence DaTA
constants and sorting of exciting lines.

program completed on 3/91 by richard a. waldo

Keyword arguments:
A -- apple
nels,zx,line,ax,ec,symb,e0,
     &exciter,fmac,phipar,caller,kn
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
from trnsprob import trnsprob


def charfdta(nels, zx, line, ax, ec, symb, e0, exciter, fmac, phipar,
             caller, kn)
      real zx(15),ax(15),fmac(15,200),exciter(200,3),ec(15),e0(15)
      integer line(15),limits(3,2)
      character*3 lines(12)
      character*1 phipar,caller
      character*7 symb(15),symbol
      logical ftest

    lines = ['Kb ', 'Ka ', 'Lg2', 'Lb3', 'Lb4', 'Lg1', 'Lb1', 'Lb2', 'La1',
               'Mg ', 'Mb ', 'Ma ']
    limits = [1, 3, 10, 2, 9, 12]

    if caller == 'C':
        nels0 = kn
        nels1 = kn
    else:
        nels0 = 1
        nels1= nels
    l = 0
    for na in range(nels0, nels1):
      # do 6000 na=nels0,nels1
        if zx[na] < 10.0:
            exciter[185 + na, 1] = l
            continue
        la = line[na]
        lea = nedge[la]
        fc5 = trnsprob(zx[na], line[na])
        for nexcit in range(1,nels):  # to 5000
            if zx[na] = zx[nexcit]:
                continue
            for k1 in range(0,2): # to 4020
                sum1=0.
                sum2=0.
                sum3=0.
                sum4=0.
                for lexcit in range(limits[k1], limits[k1 + 3]):
                    # do 4000 lexcit=limits(k1,1),limits(k1,2)
                    if ftest(e0[na], ec[na], zx[nexcit],
                             lexcit, lines[lexcit]):
                        leb = nedge[lexcit]
                        symbol = symb[nexcit][1:2], lines[lexcit], '  '
                        fc2 = effyld(zx[nexcit], leb, e0[na], 'electrons')
                        fc4 = trnsprob(zx[nexcit], lexcit)
                        fc7 = znl[leb]
                        fc8 = qe0(edge(int(zx[nexcit]), symbol[3:5]), e0[na],
                                  lines[lexcit], zx[nexcit], phipar)
                        fc = fc2*fc4*fc7*fc8
                        sum1 = sum1 + fc
                        fc1 = absionrt(zx[na], lea, xray(zx(nexcit), lexcit))
                        fc3 = effyld(zx[na], lea, xray(zx[nexcit], lexcit),
                                     'xrays    ')
                        fc6 = 1.0/ax[nexcit]
                        call macstd(zx(nexcit),zx(na),ax(na),lexcit,
            &          symb(nexcit)(1:2),symb(na)(1:2),fc9,'N','F')
                        sum2=sum2+fc1*fc2*fc3*fc4*fc5*fc6*fc7*fc8*fc9
                        sum3=sum3+edge(int(zx(nexcit)),symbol(3:5))*fc
                        xtemp=xray(zx(nexcit),lexcit)
                        sum4=sum4+xtemp*fc
                        jj=lexcit
4000        continue
            if (sum1.ne.0.) then
              l=l+1
              xedge=sum3/sum1
              xenergy=sum4/sum1
              do 4010 k2=1,nels
                call mac(xmu,int(zx(k2)),ax(k2),xenergy,iflag)
                fmac(k2,l)=xmu
4010          continue
              exciter(l,1)=sum2
              exciter(l,2)=xedge
              exciter(l,3)=(nexcit-1)*12+jj
            endif
4020      continue
5000    continue
5999    exciter(185+na,1)=l
6000  continue
      exciter(185,1)=0.
      return
      end

if __name__ == '__main__':
    pass
