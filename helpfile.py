#!/usr/bin/env python
"""-------------------------------HELPFILE-----------------------------------

author r.a. waldo 7/91
Keyword arguments:
ifile -- which help file (1-4)
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys


def helpfile(ifile):
    f = 'help%d.dat' % ifile
    with open(f, 'r')as data_file:
        data = data_file.read()
    print data


if __name__ == '__main__':
    print "File 1"
    helpfile(1)
    print "File 2"
    helpfile(2)
    print "File 3"
    helpfile(3)
    print "File 4"
    helpfile(4)
