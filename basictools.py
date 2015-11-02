#!/usr/bin/env python
"""This module contains some basic tools for input, output and usful vars.
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import numpy as np

at_els = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na',
          'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',
          'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge',
          'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo',
          'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
          'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm',
          'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
          'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb',
          'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U')
at_lines = ('Ka', 'Kb', 'Lg2', 'Lb3', 'Lb4', 'Lg1',
            'Lb1', 'Lb2', 'La1', 'Mg', 'Mb', 'Ma')


def yes_no(prompt=None, default=True):
    """return True or False based on user input
    default is returned if user input is ''
    True: y, Y, Yes, yes, YES, True, true, TRUE, t, T, 1, on, On, ON
    False: n, N, No, no, NO, False, false, FALSE, f, F, 0, off, Off, OFF
    """

    if default is None:
        prompt += "(y/n):"
    elif default is True:
        prompt += "(Y/n):"
    elif default is False:
        prompt += "(y/N):"

    true = ['Y', 'YES', 'TRUE', 'T', '1', 'ON']
    false = ['N', 'NO', 'FALSE', 'F', '0', 'OFF']
    while True:
        choice = raw_input(prompt).upper()
        if default is not None and choice == '':
            return default
        elif choice in true:
            return True
        elif choice in false:
            return False
        else:
            print "Please respond 'y' or 'n':"


def get_options(message, options=('Y', 'N'), default=None, count=False,
                casesen='Capitalize'):

    """ options = 'lines', 'els', or any tuple of strings"""
    if options == 'lines':
        options = at_lines
    if options == 'els':
        options = at_els

    stop = False
    while True:
        if casesen == 'Capitalize':
            var = raw_input(message).capitalize()
        elif casesen == 'Upper':
            var = raw_input(message).upper
        elif casesen == 'Lower':
            var = raw_input(message).lower()
        elif casesen == 'Yes':
            var = raw_input(message)

        if default is not None and var == '':
            var = default

        num = 0
        for opt in options:
            num += 1
            if var == opt:
                stop = True
                break
        if stop:
            break
    if count is True:
        return num, var
    return var


def get_nums(message, maxi=10e10, mini=0, **kwargs):
    """Get numbers from user and make sure it is in range
    Keyword Arguments:
    default -- value if just press enter i.e ''
    zeroval -- replace 0 with this value
    """

    while True:
        out = raw_input(message)
        if 'default' in kwargs and out == '':
            out = kwargs['default']
        try:
            out = float(out)
            if out > maxi or out < mini:
                print 'Must be between %d and %d.' % (mini, maxi)
            else:
                break
        except ValueError:
            print "Oops!  That was not a valid number.  Try again..."
    if 'zeroval' in kwargs and out == 0:
        out = kwargs['zeroval']
    return out


def get_string(message):
    return raw_input(message)


def out_data(header, data, width=10):
    """print lists(matrix) as table with header of column width width
    Formatted this way to make gui creation easier
    """

    if type(header) != list and type(header) != tuple:
        raise ValueError('header must be of type list or tuple')
    row_format = "{:^{width}}" * (len(header))
    print row_format.format(*header, width=width)

    if type(data) != list and type(data) != tuple:
        raise ValueError('data must be of type list or tuple')
    if type(data[0]) == list or type(data[0]) == tuple:
        for row in data:
            row_format = "{:^{width}}" * (len(row))
            print row_format.format(*row, width=width)
    else:
            row_format = "{:^{width}}" * (len(data))
            print row_format.format(*data, width=width)


def user_alert(message):
    print message


def wtfract(els, form):
    """Convert atomic formula to weight percent"""
    cnc = []
    form = np.asarray(form)
    print form
    masses = np.asarray([el.mass for el in els])
    print masses
    sum_mass = np.sum(masses*form)
    print sum_mass
    cnc = (masses*form)/sum_mass
    return cnc


if __name__ == '__main__':
    # print yes_no("Q?")
    # print get_nums("test", default=1, zeroval=-1)
    # print get_options('el', 'els')
    # out_data(('#', '#'), ((1, 1), (2, 2)))
    # out_data(('#', '#'), (1, 1), width=15)
    from atomic_element import AtomicElement as AtEl
    Si = AtEl('Si', 'Ka')
    o = AtEl('O', 'Ka')
    els = [Si, o]
    form = [1, 2]
    print wtfract(els, form)
    print 28.0/(28+16*2)
