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
from atomic_element import AtomicElement


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
                    print ('Element   Atomic Formula  X-Ray Line')
                    for j in range(1, len(data), 3):
                        print ('%s\t  %s\t  %s' % (data[j], 
                                                   round(float(data[j + 2]),
                                                               5), 
                                                   data[j + 1]))
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
            fr = []
            count = 0
            while True:
                count += 1
                print "Element %d in %s" % (count, self.name)
                name = get_options('Element Symbol:', 'els')
                line = get_options('Lowest Excited X-Ray Line:', 'lines')
                self.els.append(AtomicElement(name, line))
                if how.upper() == 'W':
                    self.els[-1].wtfract = get_nums('Weight Fraction:', 1)
                else:
                    fr.append(get_nums('Atomic Formula:', 100))
                print ("Entered so far: %s"
                       % ' '.join([o.name for o in self.els]))
                if not yes_no("Add another element? (Y/N)"):
                    break
            if how.upper() == 'A':
                wtfracs = wtfract([o.name for o in self.els], fr, [count])
                for x, _ in enumerate(wtfracs):
                    self.els[x].wtfrac = wtfracs[x]
                print ('%s compound standard: els, weight fractions:'
                       % self.name)
            print "\t".join([self.els[x].name for x
                             in range(0, len(self.els))])
            print "\t".join([str(round(o.wtfrac, 5)) for o in self.els])

            if yes_no("Is this correct?\n(restarts if no)(Y/N):", False):
                break
        data = ""
        for i in range(0, len(self.els)):
            data += "\t%s\t%s\t%s" % (self.els[i].name,
                                  self.els[i].line,
                                  self.els[i].wtfrac)
        towrite = "%s%s\n" % (self.name, data)
        with open('empastds.txt', 'a') as myfile:
            myfile.write(towrite)

if __name__ == '__main__':
    std = EpmaStd()
    print std.name
