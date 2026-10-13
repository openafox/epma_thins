#!/usr/bin/env python
# #-*-coding:utf-8-*-
"""This is my doc string.

Keyword arguments:
A -- apple
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
import numpy as np

def test(file1):
    """c string
    """
    print 'run'
    with open(file1) as fp:
        print "opened ", fp
        for i, line in enumerate(fp):
            'looping'
            if "\xe2" in line:
                print i, repr(line)
    print "done"


if __name__ =='__main__':
    file1 = "/Users/austinfox/Google_Drive/_Materials_Engr/_Project_Proposals/EPMA_thins/epma_thins/master/pap.py"
    test(file1)
