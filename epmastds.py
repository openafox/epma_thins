#!/usr/bin/env python
""" EMPA STDS

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import readline  # nice line editor for inputs
from basictools import yes_no, get_options, get_nums
from get_data import lookup
from wtfract import wtfract
import atomic_element
from atomic_element import AtomicElement


class ElementStd(atomic_element.AtomicElement):
    """ AtomicElement + wtfrac"""

    def __init__(self, name=None, line=None, atform=None):

        if name is None and line is None and atform is None:
            self.z, self.name = get_options('Please enter Element Symbol:',
                                            'els', count=True)
            self.line = get_options('Please enter Element X-Ray Line '
                                    '(Ka, Lg2, etc.):', 'lines')
            self.atform = get_nums("Weight faraction?:", 1)

        elif name is None or line is None or atform is None:
            # warn programer if improper usage
            raise ValueError('All or none of "name", "line", and "wtfrac"'
                             'must be specified.')
        else:
            self.name = name
            self.z = at_els.index(name.capitalize()) + 1
            self.line = line
            check = at_lines.index(line)
            self.atform = atform

        self.setup_vars()


class EpmaStd(object):
    """Containers for standard materials and their properties"""

    def __init__(self):
        self.name = raw_input('Enter standard name:')
        self.els = []
        if not self.load_saved():
            self.new_std()

    def load_saved(self):
        # see if exist in std doc
        with open('epmastds.txt', 'r') as data_file:
            # stdnm \t nel \t el \t wtfr \t el \t wtfr....
            for row in data_file:
                data = row.split("\t")
                if self.name.upper() == data[0].upper():
                    # Confirm std selection
                    for j in range(1, len(data), 3):
                        print ('Element   Atomic Formula\n'
                               '%s\t\t%s' % (data[j], data[j + 1]))
                    reply = yes_no('This standard ?:')
                    if reply is True:
                        for j in range(1, len(data), 3):
                            self.els.append(AtomicElement(data[j],  # name
                                                          data[j + 1],  # ln
                                                          ))
                            self.els[-1].wtfract = data[j + 2]
                        return True
                    else:
                        self.name = raw_input('Enter standard name:')
                        self.load_saved()
            # if no matches
            return False

    def new_std(self):
        how = get_options('Input %s standard composition in wt. %% (w) '
                          'or atomic formula(a) (def=a):' % self.name,
                          ('W', 'A'), 'A')
        while True:
            self.els = []
            count = 0
            while True:
                count += 1
                print "Element %d in %s" % (count, self.name)
                name = get_options('Element Symbol:', 'els')
                line = get_options('Lowest Excited X-Ray Line:', 'line')
                self.els.append(AtomicElement(name, line))
                if how.upper() == 'W':
                    self.els[-1].wtfract = get_nums('Weight Fraction:', 1)
                else:
                    fr.append(get_nums('Atomic Formula:', 100))
                print "Entered so far: %s" % ' '.join(self.els[:].name)
                if not yes_no("Add another element? (Y/N)"):
                    break
            if how.upper() == 'A':
                self.els[:].wtfrac = wtfract(self.els, fr, [count])
            print ('%s compound standard: els, weight fractions:'
                   % self.name)
            print "\t".join(self.els[:].name)
            print "\t".join([str(round(x, 5)) for x in self.els[:].wtfrac])

            if yes_no("Is this correct?\n(restarts if no)(Y/N):", False):
                break
        data = ""
        for i in range(0, len(self.els)):
            data += "\t%s\t%s" % (self.els[i].name,
                                  self.els[i].line,
                                  self.els[i].wtfrac)
        towrite = "%s\t%s\n" % (self.name, data)
        with open('standard.txt', 'a') as myfile:
            myfile.write(towrite)

if __name__ == '__main__':
    std = EpmaStd()
    print std.name
