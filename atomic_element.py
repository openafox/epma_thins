#!/usr/bin/env python
"""Atomic Element Class
container for elements and their properties
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

from basictools import (get_options, get_nums, at_els, user_alert,
                        get_data, all_or_none)


class AtomicElement(object):
    """container for elements and their properties"""
    #  #######This whole voltage (e0) thing is still very confusing
    # ######## may require some major rework depending on how it needs to work
    def __init__(self, name=None, line=None, volt=None, opt=None):
        if all_or_none(name, line, opt):
            name, line, z, volt, opt= self.user_input(volt)
            # volt can be specified if the others are not

        # Atomic Symbol
        self.name = name
        # Atomic Number
        self.z = at_els.index(name.capitalize()) + 1
        # X-Ray emmision line
        self.line = line
        # Accelerating voltage used for this element
        self.volt = volt

        self.setup_vars()

        # Setup Unknowns/initial values (Place holders)

        # Working Concentration (Maybe Kratio and weight fraction)
        self.ci = 0
        # Log previous c values (clog[0]-values, clog[1]-setby)
        self.clog = [[0, 'init']]

    def user_input(self, volt):

        z, name = get_options('Please enter Element Symbol:',
                              'els', count=True)
        # make sure the selected line is valid
        while True:
            line = get_options('Please enter Element X-Ray Line '
                               '(Ka, Lg2, etc.):', 'lines')
            xray = get_data(z, line)
            if xray > 0:
                break
            user_alert('Invalid line for this element; try again:')
        opt = self.get_opt()
        if volt is None:
            volt = get_nums('Accelerating voltage for this element?:',
                            50, 0)
        return name, line, z, volt, opt

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

        # Get element data from tables or calculation
        # Electron shell e.g. K, L1, M1, N1
        self.shell = self.get_shell()
        # Fluorescence Yields
        self.omega = self.get_omega()
        # Number of electrons in the ionized shell.
        self.znl = self.get_znl()
        # JUMP Ratios
        self.rjump = self.get_rjump()
        # EFFective Fluorescence YieLD
        self.effyld = self.get_effyld()
        # TRaNSition PROBability of line after shell ionization
        self.trnsprob = self.get_trnsprob()

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
                    if xray > 0.0:
                        self.xray = xray
                    else:
                        50000.0
                        user_alert('X-Ray emmision line not valid for'
                                   'element\n X-Ray energy set to 50MeV.'
                                   '\nPlease select another line to prevent'
                                   'anomalous calculations!')

    def get_rjump(self, z=None, shell=None):
        if all_or_none(z, shell):
            z = self.z
            shell = self.shell

        rjump = {'K': get_data(z, 'rjump1'),
                 'L1': 1.17,
                 'L2': 1.39,
                 'L3': get_data(z, 'rjump4'),
                 'M1': 1.16,
                 'M2': 1.207,
                 'M3': 1.158,
                 'M4': 1.895,
                 'M5': 1.895,
                 'N': 2.0}
        return rjump[shell]

    def get_omega(self, z=None, shell=None):
        """Fluorescence Yields"""
        if all_or_none(z, shell):
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

    def get_effyld(self, el=None, exciter='e'):
        """EFFective Fluorescence YieLD for all lines.
        3/91 r.a. waldo
        NEED TO FIND REFS FOR THIS
        """
        # need to check this xray may actually be accelerating voltage...?
        # var was xenergy in RAW code
        if el is None:
            el = self
        # Coster-Kronig Coefficients
        ck12 = get_data(el.z, 'c12')
        ck13 = get_data(el.z, 'c13')
        ck23 = get_data(el.z, 'c23')

        if el.shell[0] != 'L':
            effyld = el.omega
        elif el.shell[0] == 'L':
            rj1 = self.get_rjump(el.z, 2)
            rj2 = self.get_rjump(el.z, 3)
            rj3 = self.get_rjump(el.z, 4)
            l1 = self.get_data(el.z, 'Lb3')
            l2 = self.get_data(el.z, 'Lb1')
            l3 = self.get_data(el.z, 'La1')
            if el.shell == 'L1':
                effyld = el.omega
            elif el.shell == 'L2':
                if el.xray < l1 and el.xray > l2:
                    effyld = el.omega
                elif el.xray > l1:
                    el1 = el
                    el1.shell = 'L2'  # 3
                    omegaL2 = self.get_omega(el1)
                    if exciter == 'e':
                        effyld = omegaL2*(1.0 + ck12)
                    else:
                        effyld = omegaL2*(1.0 + ck12*(rj1 - 1.0) *
                                          rj2/(rj2 - 1.0))
            elif el.shell == 'L3':
                if el.xray < l2 and el.xray > l3:
                    effyld = el.omega
                elif el.xray < l1 and el.xray > l2:
                    effyld = el.omega*(1.0 + ck23*(rj2 - 1.0) *
                                       rj3/(rj3 - 1.0))
                elif el.xray > l1:
                    if exciter[0] == 'e':
                        effyld = el.omega*(1.0 + 0.5*ck13 +
                                           0.5*(1.0 + ck12)*ck23)
                else:
                    effyld = el.omega*(
                            1.0 +
                            ck13*(rj1 - 1.0)*rj2*rj3/(rj3 - 1.0) +
                            ck23*(rj2 - 1.0)*rj3/(rj3 - 1.0) +
                            ck12*(rj1 - 1.0)*rj2*rj3/(rj3 - 1.0))
        elif el.shell[0] == 'N':
            effyld = 0.0

        return effyld

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

    def get_trnsprob(self, z=None, line=None):
        """TRaNSition PROBability of line after shell ionization."""
        if z is None:
            z = self.z
        if line is None:
            line = self.line

        if line[0] == 'K':
            if z < 11.0:
                trnsprob = 0.01
            else:
                trnsprob = 8.0e-3*(z - 10) - 8.75e-5*(z - 10)**2
            if line == 'Ka':
                trnsprob = 1.0 - trnsprob
        elif line == 'Lg2':
            trnsprob = 0.25
        elif line == 'Lb3':
            trnsprob = 0.45
        elif line == 'Lb4':
            trnsprob = 0.3
        elif line == 'Lg1' or line == 'Lb1':
            trnsprob = (z - 60.0)*0.001 + 0.14
            if line == 'Lb1':
                trnsprob = 1.0 - trnsprob
        elif line == 'Lb2' or line == 'La1':
            trnsprob = (z - 60.0)*0.001 + 0.16
            if line == 'La1':
                trnsprob = 1.0 - trnsprob
        elif line[0] == 'M' or line[0] == 'N':
            trnsprob = 1.0

        return trnsprob

if __name__ == '__main__':
    name = 'Mg'
    line = 'Ka'
    volt = 15
    el = AtomicElement(name, line, volt)
    print 'name:', el.name
    print 'line:', el.name
    print 'mass:', el.mass
    print 'absorption edge[kV]', el.edge
    print 'electron shell', el.shell
    print 'fluorescence yields', el.omega
    print 'effective yield', el.effyld
    print 'number of electrons in the ionized shell', el.znl
    # print 'coster-kronig coefficients', el.ck
    # commented out for now, need to know how to choose which one
    print 'x-ray emmision line energy [KeV]', el.xray
    el = AtomicElement(name)
