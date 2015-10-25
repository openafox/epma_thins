#!/usr/bin/env python
"""This module contains some basic tools for input and usfull vars.
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.


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
                casesen = 'Capitalize'):

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
    if count == True:
        return num, var
    return var


def get_nums(message, maxi=10e10, mini=0, default=None):
    while True:
        out = raw_input(message)
        if default is not None and out == '':
            out = default
        try:
            out = float(out)
            if out > maxi or out < mini:
                print 'Must be between %d and %d.' % (mini, maxi)
            else:
                break
        except ValueError:
            print "Oops!  That was no valid number.  Try again..."
    return out

if __name__ == '__main__':
    # print yes_no("Q?")
    # print get_nums("test", 10, 0)
    print get_options('el', 'els')
