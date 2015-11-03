#!/usr/bin/env python
"""This is my doc string.

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
from atomic_element import AtomicElement
from basictools import get_nums, yes_no, user_alert, get_options, wtfract


class FilmLayer(object):
    """Container for properties specific to a Film Layer
    Keyword Arguments:
    els -- List of element objests
    rho -- Layer Density"""

    def __init__(self, **kwargs):
        # Is Layer of known composition and thickness?
        self.fix = False  # all layers are default not fixed
        if 'els' in kwargs:
            # List of atomic element objects in the film
            self.els = kwargs['els']
        else:
            # List of atomic element objects in the film
            self.els = self.get_els()

        if 'rho' in kwargs:
            # Layer Density
            self.rho = kwargs['rho']
        else:
            # Layer Density
            self.rho = self.get_dens()
        # Layer thickness [Angs] - unknows start as 0
        self.thick = 0.0

    def get_els(self, els=[]):
        """Get new elements to add to the layer: Returns AtomicElement object.
        """
        while True:
            els.append(AtomicElement())
            els[-1].opt = AtomicElement.get_opt
            reply = yes_no('Add another element?(Y/N):')
            if reply is False:
                break
        return els

    def get_dens(self):
        """Gets estimated layer density from user.
        Could make this a calc from xtal struct in future...
        """
        user_alert('Layer Density:\nUsed solely to convert thicknesses in '
                   'ug/cm^2 to Angstroms; has no effect on the results.')
        rho = get_nums('Approximate density (g/cm^3):', 25, 0)
        return rho

    def fixlayer(self):
        """Get fixed values for  thickness and composition"""
        # Get Thickness (layers not sub - sub has "infinite" thickness)
        if not self.fix:
            self.thick = 0.0
            return False
        self.units = get_options('+Enter fixed thicknesses in Angstroms (a)'
                                 'or ug/cm**2 (u), (def=a):', ('A', 'U'), 'A')
        if self.units == 'A':
            label = 'Angstroms'
        if self.units == 'U':
            label = 'ug/cm**2'
        self.thick = get_nums('Fixed thickness [%s]:' % label,
                              10e5, 0, 0)
        self.thick = self.thick/1.0e6   # g/cm**2
        if self.units == 'A':
                self.thick = self.thick*self.rho/100.  # cm

        # Get Composition
        if len(self.els) == 1:
            # fraction in layer is one if only 1 element
            self.els[0] = 1.0
        else:
            reply = get_options('Enter the fixed composition %s by '
                                'weight fractions (w) or by atomic '
                                'formula (def=a) :', ('W', 'A'), 'A')

            for i in range(0, len(self.els)):
                if reply == 'W':
                    for j in range(0, self.els):
                        self.els[j].wtfrac = get_nums(
                                'Weight fraction concentration'
                                'of element %s:' % self.els[j], 1.0, 0.0,
                                default=1, zeroval=1.0e-9)
                elif reply == 'A':
                        self.els[j].atform = get_nums(
                                'Atomic Formula of element %s '
                                '(i.e 1 for Si in SiO2):' % self.els[j], 100)

            if reply == 'A':
                wtfracs = wtfract([el.name for el in self.els],
                                  [el.atform for el in self.els])
                for x, _ in enumerate(wtfracs):
                    self.els[x].wtfrac = wtfracs[x]

if __name__ == '__main__':
    from atomic_element import AtomicElement as AtEl
    Si = AtEl('Si', 'Ka')
    o = AtEl('O', 'Ka')
    layer = FilmLayer(els=[Si, o], rho=2.65)
    print layer.els[0].mass
    layer2 = FilmLayer()
