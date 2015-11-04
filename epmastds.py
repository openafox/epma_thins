#!/usr/bin/env python
""" EMPA STDS

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

from basictools import (yes_no, get_options, get_nums, out_data, get_string,
                        wtfract)
from atomic_element import AtomicElement
from mac import mac


class EpmaStd(object):
    """Containers for standard materials and their properties"""

    def __init__(self, name=None):
        if name is None:
            self.name = get_string('Enter standard name:')
        self.els = []
        if not self.load_saved():
            self.new_std()

    def load_saved(self):
        # see if exist in std doc
        reply = True
        with open('epmastds.txt', 'r') as data_file:
            for row in data_file:
                data = row.split("\t")
                if self.name.upper() == data[0].upper():
                    # Confirm std selection
                    header = ['Element', 'Atomic_Formula', 'X-Ray_Line']
                    out = []
                    for j in range(1, len(data), 3):
                        out.append([data[j], round(float(data[j + 2]), 5),
                                    data[j + 1]])
                    out_data(header, out, width=15)
                    reply = yes_no('Is this data correct for standard %s?:'
                                   % self.name)
                    if reply is True:
                        for j in range(1, len(data), 3):
                            self.els.append(AtomicElement(data[j],  # name
                                                          data[j + 1],  # ln
                                                          ))
                            self.els[-1].wtfract = data[j + 2]
                        return True
                        # bounce out False using default reply
        if reply is False:
            reply = yes_no('Replace Saved Data?:')
            if reply is True:
                with open('epmastds.txt', 'r') as old_file,\
                     open('epmastds.txt', 'w') as new_file:
                    for row in old_file:
                        data = row.split("\t")
                        if self.name.upper() != data[0].upper():
                            new_file.write(row)
                self.new_std()
                return True
            else:
                self.__init__()
                return True

            # if no matches
            return False

    def new_std(self):
        how = get_options('Input %s standard composition in wt. %% (w) '
                          'or atomic formula(a) (def=a):' % self.name,
                          ('W', 'A'), 'A')
        while True:
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
                    self.els[-1].atform = get_nums('Atomic Formula:', 100)
                print ("Entered so far: %s"
                       % ' '.join([el.name for el in self.els]))
                if not yes_no("Add another element? (Y/N)"):
                    break
            if how.upper() == 'A':
                wtfracs = wtfract(self.els, [el.atform for el in self.els])
                for x, _ in enumerate(wtfracs):
                    self.els[x].wtfrac = wtfracs[x]
                print ('%s compound standard: els, weight fractions:'
                       % self.name)
            print "\t".join([self.els[x].name for x
                             in range(0, len(self.els))])
            print "\t".join([str(round(el.wtfrac, 5)) for el in self.els])

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

    def macstd(self, el1, el2, macchang, caller):
        """controls calculation of mass absorption coefficients
        for compound standards."""
        xmu = mac(el1, el2)
        if macchang == 'Y':
            mess = ('+MAC for %s %s in %s is %.4g\nChange to:(def.=current)'
                    % (el1.name, el1.line, el2.name, xmu))
            if caller == 'B':
                mess = mess + 'XRF, Compd. Stnd.:'
            if caller == 'C':
                mess = mess + 'Compound Standard:'
            if caller == 'P':
                mess = mess + 'Pure Element Stnd:'
            if caller == 'F':
                mess = mess + 'Addl. XRF, Sample:'
            xmu = get_nums(mess, 10e10, 0, xmu)
            print '+Changed to %.4g' % xmu
        return xmu

if __name__ == '__main__':
    std = EpmaStd()
    print std.name
