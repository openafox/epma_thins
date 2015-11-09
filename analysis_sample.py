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

    def get_layers_depth(self):
        depth = 0.0
        for i in range(len(self.layers), 0):
            depth += self.layers[i].thick
            self.layers[i].depth = depth

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
        """get take-off angle"""
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

    def get_chiov(self, el1, layindex):
        """chi OVerlayer
        finds the chi (wt. fraction averaged mass absorption coefficient)
        of an overlayer for an element 'el1' in an underlying layer (layindex)
        """

        chiovl = list(np.zeros((0, len(self.layers) - layindex)))
        # surface to current layer
        for lar, el2 in self:
            if self.layers.index(lar) > layindex:  
                # only calculate to current layer
                chiovl[i] += el2.c1 * el1.mac[el2.name]
        return chiovl

    def get_macs(self):
        """Get Mass Absorption Coefficients
        Each element has a mac for its x-ray emmision by every other element.
        """
        for lar1, el1 in self:
            el1.mac = {}
        for (lar1, el1), (lar2, el2) in itertools.product(self, repeat=2):
                # print el1.name, el2.name
                if not el2.name in el1.mac:
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

    def __iter__(self):
        self.j = 0
        self.i = len(self.layers) - 1
        return self

    def next(self):
        layindex = self.i
        elindex = self.j
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
    Si = AtEl('Si', 'Ka', 15, 's')
    o = AtEl('O', 'Ka', 15, 's')
    ti = AtEl('Ti', 'Ka', 15, 's')
    layer1 = FL(els=[Si, o], rho=2.65)
    layer2 = FL(els=[ti, o], rho=4.23)
    sample = AnalysisSample(toa=40, volts=[15], layers=[layer1, layer2],
                            phimodel='E')
    print 'defs done'
    for lay, el in sample:
        print ('Element:', el.name, ' Macs', el.mac, 'layer',
               sample.layers.index(lay))
    print ""
    for lay, el in sample:
        print 'Element:', el.name, 'Density', lay.rho
