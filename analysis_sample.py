#!/usr/bin/env python
"""This is my doc string.

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
from basictools import get_nums

class AnalysisSample(object):
    """ Container for all analysis properties"""

    def __init__(self):
        # Take Off Angle
        self.toa = self.get_toa
        # List of accelerating voltages used on all elements
        # gmrfilm defined for each
        self.volts = []
        get_nums('\tEnter Eo for element (def.=15'
                                          '):', 58, 0, 15))
        # List of film layers
        self.layers = self.get_filmlayers()
        # iteration vars
        self.i = 0
        self.j = 0


    def get_filmlayers(self)
        
    def check_volts(self):
        for el in self.layers:
            if el.ecr <= self.volts*0.9:
                #here#
                break
                print 'Overvoltage ratio for %s (%5.1f keV Ec) with Eo='
                      '%5.1f keV is too low!\nTry Again:'
                      % (els[i], ecr, e0[i])

    def get_toa():
        """get take-off angle"""
        self.toa = get_nums('Enter take-off angle (default=40 degrees):',
                            90, 0, 40)
    def fixlayers():
        
        for i in range(0, len(self.els)):
            if i > 0:
                mess = '+Fix composition and thickness of layer %d? (def=n):' % i)
            if i == 0:
                mess = '+Fix composition of substrate? (def=n):')
            ifix[i] = yes_no(mess, False)

                    if i > 0:
                        label = 'of layer %d' % i
                    else:
                        label = 'of the substrate'

    # may be useful but add to 
    def newline(nels, els, e0, lines):
        """Queries operator for a new line for an element
        because the old line was not excited at the new operating potential.

        program completed 5/91 by richard a. waldo

        Keyword Arguments:
        nels -- Number of Elements
        els -- List of element symbols
        lines -- List of element lines used
        e0 -- ??Original X-Ray energies list
        Return:
        ec -- Element Line X-Ray Energy [KeV]
        lines -- updated
        """
        for i in range(0, nels):
            (z, mass, ecr) = lookup(els[i], lines[i])
            while True:
                if ecr > e0[i]*0.9:
                    lines[i] = get_options(
                            'Overvoltage ratio for %s (%5.1f keV Ec) with Eo='
                            '%5.1f keV is too low!\nEnter another line (e.g.'
                            ' La1, Ka):' % (els[i], ecr, e0[i]), options, 'Ka')
                elif ecr <= 0:
                    lines[i] = get_options('Invalid line for this element;'
                                        'try again:', options, 'Ka')
                else:
                    break
                (z, mass, ecr) = lookup(els[i], lines[i])
            ec[i] = ecr
        return (ec, lines)


    def __iter__(self):
        return self

    def next(self):
        nlays = len(self.layer)
        lnel = len(self.layer[nlays-1].element)
        try:
            layer = self.layer[self.i]
            element = self.layer[self.i].element[self.j]
        except IndexError:
            self.i = 0
            self.j = 0
            raise StopIteration
        if self.j == len(self.layer[self.i].element) - 1:
            self.i = self.i + 1
            self.j = 0
        else:
            self.j = self.j + 1
        return layer, element
    

if __name__ == '__main__':
    pass
