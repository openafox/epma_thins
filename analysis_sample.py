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


class AnalysisSample(object):
    """ Container for all analysis properties
    Keyword Arguments(give all or none):
    toa -- Take Off Angle
    volts -- list of accelerating voltages used
    layers -- list of FilmLayer objects
        (these contain an els attribute -- list of AtomicElement objects)
    """

    def __init__(self, **kwargs):
        if len(kwargs) == 3:
            # Take Off Angle
            self.toa = kwargs['toa']
            # List of accelerating voltages used on all elements
            # gmrfilm defined for each
            self.volts = kwargs['volts']
            # List of film layers
            self.layers = kwargs['layers']
        else:
            # Take Off Angle
            self.toa = self.get_toa
            # List of accelerating voltages used on all elements
            # gmrfilm defined for each
            self.volts = self.get_volts = []
            # List of film layers
            self.layers = self.get_layers([])

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

    def __iter__(self):
        self.j = 0
        self.i = len(self.layers) - 1
        return self

    def next(self):
        layindex = self.i
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
        return layindex, layer, el


if __name__ == '__main__':
    from atomic_element import AtomicElement as AtEl
    from film_layer import FilmLayer as FL
    Si = AtEl('Si', 'Ka')
    o = AtEl('O', 'Ka')
    ti = AtEl('Ti', 'Ka')
    layer1 = FL(els=[Si, o], rho=2.65)
    layer2 = FL(els=[ti, o], rho=4.23)
    sample = AnalysisSample(toa=40, volts=[15], layers=[layer1, layer2])
    print 'o'
    for i, lay, el in sample:
        print 'Element:', el.name, ' Mass', el.mass, 'layer', i
    print ""
    for i, lay, el in sample:
        print 'Element:', el.name, 'Density', lay.rho, 'layer', i
    
