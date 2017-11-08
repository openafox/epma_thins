#!/usr/bin/env python
"""------------------------------STANDARD---------------------------------

calculation of theoretical standard intensities
for both pure element and compound standards
4/87 r.a. waldo
3/88 to include compound standards, r.a. waldo

Keyword arguments:
A -- apple
  nels,symb,ec,e0,line,fconts,fchars
     &,stds0,stds1,stds2,prims,phipar,toa,macchang,icmd8
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
from epmastd import epmastd
from atomic_element import AtomicElement
from analysis_sample import AnalysisSample
import pap
import phirzeq


    def get_macs(self):
        """Get Mass Absorption Coefficients
        Each element has a mac for its x-ray emmision by every other element.
        """
        for lar1, el1 in self:
            el1.mac = {}
        # for all combinations of elements
        for (lar1, el1), (lar2, el2) in itertools.product(self, repeat=2):
                # print el1.name, el2.name
                if el2.name not in el1.mac:
                    el1.mac[el2.name] = Mac(el1, el2)
def standard(samp, macchang):
    # what is prims?? (total els, 2)

    ff = False
    for lay_i, el_i in samp:
        c[1] = 1.0
        dumline(1)=line(i)
        ecdum(1)=ec(i)
        symbol=symb(i)
        dumsym(1)=symb(i)

        xmu = macstd(el_i, el_i, macchang,'p')
        mcstd.append(xmu)
        mcstd(1,1)=xmu
        idum(1)=1
        fc1 = el.effyld(el)  # need to check if voltage will need to be spec
        fc2 = el.trnsprob
        fc3 = 1.0/el.mass
        fc4 = el.znl
        fc5 = self.qe0(el, el.volt, self.phimodel)
        const = fc1*fc2*fc3*fc4*fc5
            if phimodel == 'e':
               (a_1, b_1, a_2, r_m, r_c, r_x, dza) = pap.pap(samp, lay_i,
                                                             el_i, 'E')
                # E - pure element standard
                el_i.prims1 = dza
            else:
                dum(1)=0.
                alpha, beta, gamma, phi0 = phirzeq.phirzeq(samp, el_p, lay_p,
                                                           ff)
                call phirzeq(idum,1,1,dum,al,b,g,p,z,a,c,e0(i),ec(i),phipar,line(i),ff,zz)
                dza = phirzeq.integral(samp, el_p, lay_p,  8, al, b, g, p, xmu)
                call integral(dum,8,al,b,g,p,xmu,dum,prims(i,1),toa)
            el_i.prims1 = const*dza
            std=stds0(i)
            fchars(i,1)=0.

            # If continuum corrections are included do it
            if (icmd8.eq.'y') then
            call contflr(1,idum,1,z,ec(i),dumline,a,
        &    e0(i),c,fconts(i,1),phipar,toa,mcstd,delta,'b')
            endif
c
c        if a compound standard
c
        if (std.gt.0) then
          l=1
          do 50 ij=1,15
            if (stds1(l,1,std).le.0.) goto 51
            if (symb(i)(1:2).eq.stds2(l,std)) k=l
            c(l)=stds1(l,1,std)
            z(l)=stds1(l,2,std)
            a(l)=stds1(l,3,std)
            l=l+1
            ecdum(l)=0.
50        continue
51        ecdum(k)=ec(i)
          do 70 m=1,l-1
            call macstd(z(k),z(m),a(m),line(i),stds2(k,std),
     &      stds2(m,std),xmu,macchang,'c')
            mcstd(k,m)=xmu
            stds(m)=stds2(m,std)
            dumsym(m)(1:2)=stds(m)
70        continue

          fchars(i,2)=0.
          idum(1)=l-1
          dumline(k)=line(i)
          ecdum(k)=ec(i)
          dumsym(k)=symb(i)
          if (phipar.eq.'e') then
            call pap(dum,k,idum,l-1,z,a,c,e0(i),ec(i),line(i),
     &       toa,mcstd,'b',za,a1,a2,b1,rc,rm,rx,symbol,'c',zz)
          else
            call phirzeq(idum,l-1,k,dum,al,b,g,p,z,a,c,
     &      e0(i),ec(i),phipar,line(i),ff,zz)
            ch=0.
            do 75 j=1,l-1
75            ch=ch+c(j)*mcstd(k,j)
            call integral(dum,8,al,b,g,p,ch,dum,za(1),toa)
          endif
          prims(i,2)=const*c(k)*za(1)
          call charfdta(idum(1),z,dumline,a,ecdum,dumsym,e0,
     &    exciter,fmac,phipar,'c',k)
          call charflr(k,idum,idum(1),'b',z,a,dumsym,
     &    e0(i),c,dum,fchars(i,2),phipar,toa,mcstd,exciter,fmac)
          if (icmd8.eq.'y') then
            call contflr(k,idum,idum(1),z,ec(i),dumline,a,
     &      e0(i),c,fconts(i,2),phipar,toa,mcstd,delta,'b')
          endif
        else
          prims(i,2)=prims(i,1)
          fchars(i,2)=fchars(i,1)
          fconts(i,2)=fconts(i,1)
        endif
100   continue
      return
      end
if __name__ == '__main__':
    pass

