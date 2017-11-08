#!/usr/bin/env python
"""------------------------------FTest------------------------------------
Fluorescence Test; test if characteristic x-ray fluorescence
can occur for the given elements and beam voltage.   6/88 r.a.waldo


Keyword arguments:
e0 -- ??
eci -- ??
zb -- ??
lb -- ??
line -- ??
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys


def ftest(e0, eci, zb, lb, line)
    xb = edge(int(zb), line)
    if e0 > xb:
        if eci < xray(zb, lb) and 2.0*eci > xray(zb,lb):
            return True
    return False
if __name__ == '__main__':
    pass

