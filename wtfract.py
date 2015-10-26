#!/usr/bin/env python
"""-----------------------------WTFRACT---------------------------------
a subroutine to convert atomic formula
to atomic percent and weight percent
4/87 raw.

Arguments:
els -- list of elements
lines -- list of lines
form -- list of formula units ex SiO2 [1, 2]
nel -- list of elements per layer
Return:
cnc -- list of weight fractions
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

from get_data import lookup
import layrelem
import numpy as np

def wtfract(els, form, nel):
    cnc = [] 
    masses = []
    form = np.asarray(form)
    j = 0 # element counter

    for k in range(0, len(els)):
        (z, mass) = lookup(els[k])
        masses.append(mass)

    masses = np.asarray(masses)
    for i in range(0, len(nel)):
        sumwp = np.sum(form[j:j + nel[i]]*masses[j:j + nel[i]])
        for k in range(0, nel[i]):
            # print els[j], masses[j]
            cnc.append(masses[j]*form[j]/sumwp)
            j = j + 1  # element counter
    print cnc
    return cnc

if __name__ == '__main__':
    els = ['Si', 'O', 'Mg', 'O']
    form = [1, 2, 1, 1]
    nel = [2, 2]

    print wtfract(els, form, nel)
    form = [2, 4, 0.5, 0.5]
    nel =[2]
    print wtfract(els, form, nel)
    print 28.0/(28+16*2)
