#!/usr/bin/env python
"""------------------------Help-----------------------

Prints help information
5/93 r.a. waldo


Keyword arguments:

"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
from helpfile import helpfile


def gmrhelp():
    hlp = ' '
#     write(*,7001)
#7001 format (' Type (H)elp for information on the program, otherwise hi
#    &t Enter key: ')
#     read 990,hlp
    if hlp == 'h' or hlp == 'H':
        helpfile(1)
     
      return
if __name__ == '__main__':
    pass

