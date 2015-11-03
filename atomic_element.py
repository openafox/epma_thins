#!/usr/bin/env python
"""Atomic Element Class
container for elements and their properties
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

from basictools import get_options, at_els, user_alert, get_data


class AtomicElement(object):
    """container for elements and their properties"""

    def __init__(self, name=None, line=None):
        if name is None and line is None:
            # Elemet Name and Line
            self.name, self.line, self.z = self.get_el_nd_ln()
        elif name is None or line is None:
            # warn programer if improper usage
            raise ValueError('Both or neither "name" and "line" must be '
                             'specified.')
        else:
            self.name = name
            self.z = at_els.index(name.capitalize()) + 1
            self.line = line
            ecr = get_data(self.z, self.line)
            if ecr < 0:
                """need to figure out how to make this ok for gui"""
                user_alert('Invalid line for this element; try again:')

        self.setup_vars()

    def get_el_nd_ln(self):

        z, name = get_options('Please enter Element Symbol:',
                              'els', count=True)
        # make sure the selected line is valid
        while True:
            line = get_options('Please enter Element X-Ray Line '
                               '(Ka, Lg2, etc.):', 'lines')
            ecr = get_data(z, line)
            if ecr > 0:
                break
            user_alert('Invalid line for this element; try again:')
        return name, line, z

    def get_opt(self):
        """Options for analysis
        s -- Determined by stoichiometry.
        c -- analyzed with a compound std.
        d -- analyzed by difference (bulk only).
        E -- analyzed with pure element standard.
        """
        return get_options('Options((S)toic, (C)omp, (D)iff, (E)lem):',
                           ['S', 'C', 'D', 'E'], 'C')

    def setup_vars(self):
        """Get element data"""

        # Electron shell e.g. K, L1, M1, N1
        self.shell = self.get_shell()
        # Fluorescence Yields
        self.omega = self.get_omega()
        # Number of electrons in the ionized shell.
        self.znl = self.get_znl()
        with open('atomic_data.txt', 'r') as data_file:
            for row in data_file:
                data = row.split("\t")
                if int(data[0]) == 0:
                    cols = data
                if self.z == int(data[0]):

                    # Element Mass
                    self.mass = float(data[cols.index('Mass')])

                    edge = float(data[cols.index(self.shell)])
                    # Absorption Edge [KeV]
                    self.edge = edge if edge > 0.0 else 50000.0

                    # Coster-Kronig Coefficients
                    # self.ck = float(
                    #        data[cols.index(get_options('ck (c12, c13, c23):',
                    #                                    ('c12', 'c13', 'c23'),
                    #                                    casesen = 'Lower'
                    #                                    ))])

                    xray = float(data[cols.index(self.line)])
                    # X-Ray Emmision Line Energy [KeV]
                    self.xray = xray if xray > 0.0 else 50000.0

                    rjump = {'K': float(data[cols.index('rjump1')]),
                             'L1': 1.17,
                             'L2': 1.39,
                             'L3': float(data[cols.index('rjump4')]),
                             'M1': 1.16,
                             'M2': 1.207,
                             'M3': 1.158,
                             'M4': 1.895,
                             'M5': 1.895,
                             'N': 2.0}
                    # JUMP Ratios
                    self.rjump = rjump[self.shell]

    def get_omega(self):
        """Fluorescence Yields"""
        z = self.z
        shell = self.shell
        omega = 0.0
        if shell == 'K':  # Bambynek etal., Rev. Mod. Physics. 1972 pg 757
            # Could updat to Bambynek 1984
            d = (0.015 + 0.0327*z - 0.64e-6*z**3)**4
            omega = d/(1.0 + d)
        elif shell == 'M4' or shell == 'M5':  # ???? REF??
            omega = 0.68e-9*(z - 13)**4
        elif shell == 'N1':  # ??? Ref??
            omega = 1.0e-9*(z - 13)**4
        elif shell == 'M2' or shell == 'M3' or shell[0] == 'N':  # 5
            omega = 0.0
        else:  # Krause, J. Phys. Chem. Ref. Data Vol. 8, No. 2, 1979, p. 307.
            if shell == 'M1':
                shell = 'L3'
            shell = shell.lower()
            omega = get_data(z, shell)
        return omega

    def get_shell(self):
        """return shell for a line"""
        shell = 'K'
        line = self.line
        if line[0] == 'K':
            shell = 'K'
        elif line[0] == 'L':
            if line[1:3] == 'g2' or line[1:3] == 'b3' or line[1:3] == 'b4':
                shell = 'L1'
            elif line[1:3] == 'g1' or line[1:3] == 'b1':
                shell = 'L2'
            elif line[1:3] == 'b2' or line[1:3] == 'a1':
                shell = 'L3'
        elif line[0] == 'M':
            if line[1] == '1' or line[1] == '2':
                shell = line
            if line[1] == 'g':
                shell = 'M3'
            if line[1] == 'b':
                shell = 'M4'
            if line[1] == 'a':
                shell = 'M5'
        elif line[0] == 'N':
            shell = line
        return shell

    def get_znl(self):
        """Number of electrons in the ionized shell"""
        znl = 0
        shell = self.shell
        if shell[0] != 'K':
            if shell[1] == 3 or shell[1] == 4:
                znl = 4
            elif shell[1] == 5:
                znl = 6
            else:
                znl = 2
        else:
                znl = 2
        return znl

if __name__ == '__main__':
    name = 'Mg'
    line = 'Ka'
    el = AtomicElement(name, line)
    print 'name:', el.name
    print 'line:', el.name
    print 'mass:', el.mass
    print 'absorption edge[kV]', el.edge
    print 'electron shell', el.shell
    print 'fluorescence yields', el.omega
    print 'number of electrons in the ionized shell', el.znl
    # print 'coster-kronig coefficients', el.ck
    # commented out for now, need to know how to choose which one
    print 'x-ray emmision line energy [KeV]', el.xray
    el = AtomicElement(name)
