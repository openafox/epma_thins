#!/usr/bin/env python
"""This is my doc string.

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
from basictools import get_nums, yes_no, get_options
from film_layer import FilmLayer
import numpy as np
from mac import Mac
import itertools


class AnalysisSample(object):
    """ Container for all analysis properties
    Keyword Arguments(give all or none):
    toa -- Take Off Angle
    volts -- list of accelerating voltages used
    layers -- list of FilmLayer objects
        (these contain an els attribute -- list of AtomicElement objects)
    """
    #need a samp.mode B or F
    def __init__(self, **kwargs):
        if len(kwargs) == 4:
            # Take Off Angle
            self.toa = kwargs['toa']
            # List of accelerating voltages used on all elements
            # gmrfilm defined for each element
            self.volts = kwargs['volts']
            # List of film layers
            self.layers = kwargs['layers']
            # ph(rz) model type
            self.phimodel = kwargs['phimodel']
        elif len(kwargs) > 0:
            raise ValueError('all or none of kargs must be specified.')
        else:
            # Take Off Angle
            self.toa = self.get_toa
            # List of accelerating voltages used on all elements
            # gmrfilm defined for each
            self.volts = self.get_volts = []
            # List of film layers
            self.layers = self.get_layers([])
            # ph(rz) model type
            self.phimodel = self.get_phimodel
        # cosecant of take off angle
        self.csctheta = 1.0/np.sin(self.toa*180/np.pi)
        # Mass Absorption Coefficients
        self.get_macs()

    def get_layers(self, layers=[]):
        """Get new layers to add to the sample:
        Returns list of FilmLayer objects.
        """
        while True:
            layers.append(FilmLayer())
            reply = yes_no('Add another layer?(Y/N):')
            if reply is False:
                break
        return layers

    def update_layers_depth(self):
        """find depth (to top) of burried layers."""
        depth = 0.0
        for i in range(len(self.layers) - 1, 0, -1):
            self.layers[i].depth = depth
            depth += self.layers[i].thick

    def get_phimodel():
        """gets desired phi(rz) model from user"""
        phimodel = get_options(
                'Choices of phi(rz) models are(default=E):\n'
                '\t(B)\tBastin\'s              Scanning (1986)\n'
                '\t(C)\tBastin\'s              Scanning (1990)\n'
                '\t(E)\tPouchou, Pichoir (PAP) Scanning (1990)\n'
                '\t(P)\tPackwood\'s            MAS      (1986):\n',
                ('B', 'C', 'E', 'P'), 'E')
        return phimodel

    def get_volts(self):
        """Get accelerating voltages"""
        volts = []
        while True:
            volts.append(get_nums('Accelerating Voltage [KeV]?:', 60, 0))
            reply = yes_no('Add another voltage?(Y/N):')
            if reply is False:
                break
        return volts

    def check_volts(self):
        for layer, el in self:
            while True:
                if el.edge <= min(self.volts)*0.9:
                    break
                print ('Overvoltage ratio for %s (%5.1f keV Ec) with Eo='
                       '%5.1f keV is too low!\nTry Again:'
                       % (el.name, el.edge, min(self.volts)))
                el.line = get_options('Lower Energy X-Ray Line for %s(%s)?'
                                      '(Ka, Lg2, etc.):' % (el.name, el.line),
                                      'lines')
                el.setup_vars()

    def get_toa(self):
        """get take-off angle from user"""
        return get_nums('Enter take-off angle (default=40 degrees):',
                        90, 0, 40)

    def fixlayers(self):
        for i in range(0, len(self.layers)):
            if i > 0:
                mess = ('Fix composition and thickness of layer %d? (def=n):'
                        % i)
            if i == 0:
                mess = 'Fix composition of substrate? (def=n):'

            self.layers[i].fix = yes_no(mess, False)
            self.layers[i].fixlayer()

    def find_ovl_macs(self, lar_p, el_p):
        """finds the wt. fraction averaged mass absorption coefficient of
        overlayers for an element 'el_p' in an underlying layer 'lar_p'
        """

        layindex = self.layers.index(lar_p)
        el_p.ovl_macs = [0 for i in self.layers]  # create list of zeros
        # surface to current layer
        for lar_i, el_i in self:
            i = self.layers.index(lar_i)
            if i >= layindex:
                # only calculate to current layer leave rest as zero
                el_p.ovl_macs[i] += el_i.c1 * el_p.mac[el_i.name]

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

    def qe0(self, el, volt, phimodel=None):
        """Ionization cross section at operating potential
        3/91 r.a. waldo
        Modified 10/91 according to Pouchou and Pichoir
        in "Electron Probe Quantitation" Plenum Press (1988 NBS Workshop)
        """
        if phimodel is None:
            phimodel = self.phimodel
        z = el.z
        line = el.line
        u0 = volt/el.edge
        if phimodel == 'B' or phimodel == 'P':
                mparam = 0.8
        else:
            if line[0] == 'L':
                mparam = 0.82
            if line[0] == 'M':
                mparam = 0.78
            if line[0] == 'K':
                if z > 30:
                    mparam = 0.86
                else:
                    mparam = 0.86 + 0.12*np.exp(-1.0*(z/5.0)**2)
            # 39229.=pi x e^4 ; e is the electron charge
            # b=.76
        q = np.log(u0)/el.edge**2/u0**mparam
        if phimodel == 'B' or phimodel == 'P':
                qe0 = 39229.0*0.76*q
        elif phimodel == 'C' or phimodel == 'E':
            if line[0] == 'K':
                qe0 = 1e-20*3.8*mparam*6.023e23*q
            if line[0] == 'L':
                qe0 = 1e-20*5.7*mparam*6.023e23*q
        if line[0] == 'M':
            qe0 = 39229.0*mparam*q
        return qe0

    def get_valences(self):
        """input of valences when one element is analyzed by stoichiometry
        May want to make this a selection in gui version...
        """
        for lar, el in self:
            if el.opt == 'S':
                lar.stoich = el
            else:
                lar.stoich = 0
        for lar, el in self:
            i = self.layers.index(lar)
            if lar.fix:
                continue
            if lar.stoich == 0:
                continue
            if el.valence == -1:
                el.valence = get_nums('Enter valence for layer %d element %s:'
                                  % (i, el.name), 20, 0)

    def calc_stoich(self):
        """calculates weight and atomic fractions of all elements given
        valences if one element is analyzed by stoichiometry
        """
        self.get_valences()
        for lar in self.layers:
            if lar.fix:
                continue
            if lar.stoich != 0:
                el_s = lar.stoich
                stc = 0.
                for el in lar.els:
                    if lar.stoich != el:
                        el.atpc = el.c1 / el.mass
                        stc += (el.atpc * el.valence / el_s.valence)
                el_s.atpc = abs(stc)
                el_s.c1 = el_s.atpc * el_s.mass
                tsum = np.sum(np.asarray([x.c1 for x in lar.els]))
                for el in lar.els:
                    el.c3 = el.c1
                    el.c1 = el.c1/tsum  # normalizing concentrations??

    def calc_thick0(self):
        """Calculates starting film thicknesses for the iteration procedure.
        Based on RAW's MAS 88 paper and the strtthick code in GMRFilm
        """
        self.calc_stoich()

        for lar, el in self:
            li = self.layers.index(lar)
            r = 6.5e-6 * el.volt**1.70  # from heinrich p419
            # (Castaing, R., Adv. Elec. Phys. 13, 317 (1960)
            # why not use 7.0 ... 1.65??   from MAS 88 paper eq 16 (Heinrich
            # derivation p419) not sure why to chose on over the other, yet
            if lar.fix or lar.stoich == el.name:
                continue

            factor = 0.0
            if li > len(self.layers):
                # if burried layer
                self.find_ovl_macs(el, li)
                for j in range(len(self.layers), li, -1):
                    factor += self.layers[j].thick * el.ovl_macs[j]

            self.update_layers_depth()
            #changred kr to c1  need to check
            if self.layers.index(lar) == len(self.layers):
                # if surface layer
                lar.thick += r*(1.0 - np.sqrt(1.0 - el.c1))
                # from MAS 88 paper eq 16 (sorta, wrong in paper)
                # it is correct however and is derived in docs
                chi = 1.0  # Abs Coeff of over layers
            elif self.layers.index(lar) == 0:
                # if substrate
                chi = np.exp(-self.csctheta * factor)  # Abs Coeff overlayers
            else:
                # if burried layer
                t1 = np.sqrt(1.0 - el.c1 - 2.0*lar.depth/r + (lar.depth/r)**2)
                lar.thick += (r*(1.0 - t1) - lar.depth)
                # from MAS 88 paper eq 16 (sorta, wrong in paper)
                # it is correct however and is derived in docs
                chi = np.exp(-self.csctheta * factor)  # Abs Coeff overlayers
            el.c1 = el.c1/chi
            el.clog.append([el.c1, 'thick0'])

    def __iter__(self):
        """itterate from top layer down"""
        self.j = 0
        self.i = len(self.layers) - 1
        return self

    def next(self):
        try:
            layer = self.layers[self.i]
            el = self.layers[self.i].els[self.j]
        except IndexError:
            raise StopIteration
        if self.i < 0:
            raise StopIteration
        if self.j == len(self.layers[self.i].els) - 1:
            self.i = self.i - 1
            self.j = 0
        else:
            self.j = self.j + 1
        return layer, el


if __name__ == '__main__':
    from atomic_element import AtomicElement as AtEl
    from film_layer import FilmLayer as FL
    el1 = AtEl('Si', 'Ka', 15, 'E')
    el2 = AtEl('O', 'Ka', 15, 'S')
    el3 = AtEl('Ti', 'Ka', 15, 'E')
    el4 = AtEl('O', 'Ka', 15, 'E')
    layer1 = FL(els=[el1, el2], rho=2.65)
    layer2 = FL(els=[el3, el4], rho=4.23)
    samp = AnalysisSample(toa=40, volts=[15], layers=[layer1, layer2],
                            phimodel='E')
    print 'defs done'
    for lay, el in samp:
        print ('Element:', el.name, ' Macs', el.mac, 'layer',
               samp.layers.index(lay))
    print ""
    for lay, el in samp:
        print 'Layer:', samp.layers.index(lay), 'thickness:', lay.thick
        print 'Element:', el.name, 'c1', el.c1
    it = 0
    while True:
        it += 1
        thick1 = samp.layers[1].thick
        samp.calc_thick0()
        for lay, el in samp:
            print 'Layer:', samp.layers.index(lay), 'thickness:', lay.thick
            print 'Element:', el.name, 'c1', el.c1
        print ""
        if it > 20:
            print "done it"
            break
        if abs(thick1 - samp.layers[1].thick) < 0.0001:
            print "done diff"
            break
