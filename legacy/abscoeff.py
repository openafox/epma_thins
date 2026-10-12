#!/usr/bin/env python
""" -------------------------------ABSCOEFF---------------------------------

This program controls the calculation (via subroutines mu and mac)
of mass absorption coefficients.  The option is present for operator
input of the mac's.  Used for layer and substrate elements in a
thin film system.
author r.a. waldo 4/87
changed to Heinrich's ICXOM-11
Mass Absorption Coefficient Equations 10/88

Keyword arguments:
nel -- ??
nels -- ??
mac -- Mass Absorption Coeficient
symb -- ??
zx -- ??
ax -- ??
line -- ??
macchange -- ??
m -- ??
Return:
?? -- ??
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import mu
import layrelem
import layr


def abscoeff(nel, nels, mac, symb, zx, ax, line, macchange, m):
    #      real mac(15,15),zx(15),ax(15)
    # defining matrix mac 15x15  zx array 15 long
    #      integer line (15),nel(7)
    #      character*7 symb(15),tempabs1*12,macchange*1,m*1

    for it1 in range(1, nels):
        (nea, la) = layrelem.layrelem(nel, i1)
        #     if (macchange.ne.'Y') goto 30
        if macchange == 'Y':
            if m == 'B':
                print '+Mass absorption coefficients (def=no change)  :'
            else:
                for kk in range(1, 6):
                    if la == kk and nea == 1:
                        print '+Mass absorption coefficients for %6s %1d '\
                            'elements (def=no change)  :' % ('layer ', kk)
                        print '  '
    # 20 ####
                if la == 7 and nea == 1:
                    print '+Mass absorption coefficients for substrate '\
                            'elements (def=no change):'
                    print '  '
    # 30
        for it2 in range(1, nels):
            lb = layr.layr(nel, it2)
            if la < lb:
                continue
            (mu, zx[i1], zx[i2], line[i1], ax[i2]) = mu.mu(mu, zx[i1], zx[i2],
                                                           line[i1], ax[i2])
            if not (m == 'B' and i1 == i2) or macchange != 'Y':
                if m == 'B':
                    print '+MAC for %s in %s is : %9.2f' % (symb[i1][1:5],
                                                            symb[i2][1:2], xmu)
                elif la != 7 and lb != 7:
                    print '+MAC for layer %d  %s in layer %d  %s is : %9.2f'\
                        % (la, symb[i1][1:5], lb, symb[i2][1:2], xmu)
                elif la == 7 and lb != 7:
                    print '+MAC for substrate %s in layer %d  %s is : %9.2f'\
                        % (symb[i1][1:5], lb, symb[i2][1:2], xmu)
                elif la == 7 and lb == 7:
                    print '+MAC for substrate %s in substrate %s is : %9.2f'\
                        % (symb[i1][1:5], symb[i2][1:2], xmu)
                tempabs1 = raw_input('Change to? : ')
                if tempabs1 != ' ':
                    for i in range(0, 11):
                        if tempabs1[i] == '.':
                            break
                        if tempabs1[i] == ' ':
                            tempabs1[i] = '.'
                            break
                # 40
                # 50
                # print '+Changed to :%10.2f'9002,xmu
            mac[i1, i2] = xmu
        #  150    continue
    #  200    continue
    return


if __name__ == '__main__':
    pass
