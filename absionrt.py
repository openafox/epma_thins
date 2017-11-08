#!/usr/bin/env python
"""------------------------------absionrt--------------------------------
RaTio of IONizations to total ABSorptions
3/91 r.a. waldo

Keyword arguments:
z -- Atomic Number or Element Symbol
jedge -- ??
xenergy -- ??
Return:
absionrt -- Ionizations/Absorptions
"""
# Copyright 2015 Austin Fox
# Program is distributed under the terms of the
# GNU General Public License see ./License for more information.

import sys
from rjump import rjump
from edge import edge


def absionrt(z, jedge, xenergy):
    # goto (1,2,2,2,1,1,1,1,1,3,3,3)jedge
    #       1,2,3,4,5,6,7,8,9,10
    j = jedge
    if j == 1 or (j >= 5 and j <= 9):
        absionrt = (rjump(z, jedge) - 1.0)/rjump(z, jedge)
        return absionrt
    if j >= 2 and j <= 4:
        jr1 = rjump(z, 2)
        jr2 = rjump(z, 3)
        jr3 = rjump(z, 4)
        l1 = edge(z, 'Lb3')
        l2 = edge(z, 'Lb1')
        l3 = edge(z, 'La1')
        if jedge == 2:
            absionrt = (jr1 - 1.0)/jr1
        elif jedge == 3:
            if xenergy < l1 and xenergy > l2:
                absionrt = (jr2 - 1.0)/jr2
            elif xenergy > l1:
                absionrt = (jr2 - 1.0)/jr1/jr2
        elif jedge == 4:
            if xenergy < l2 and xenergy > l3:
                absionrt = (jr3 - 1.0)/jr3
            elif xenergy < l1 and xenergy > l2:
                absionrt = (jr3 - 1.0)/jr2/jr3
            elif xenergy > l1:
                absionrt = (jr3 - 1.0)/jr1/jr2/jr3
        return absionrt

    if j > 9:
        absionrt = 1.0
        return absionrt

if __name__ == '__main__':
    print(absionrt(12, 2, 3))
