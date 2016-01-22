# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 12:05:02 2014

@author: derek
"""
import io
import decomposition as dec
import sys
import domain


if __name__ == '__main__':
    """
    example: python analyze_hgb.py <hgb-filename>
    """
    if len(sys.argv) == 1:
        exit(0)
    print sys.argv[1]

    d = io.readForCBin(sys.argv[1])
    print d

    d.disp()

    io.writeHGA(d, "%s.hga" % (sys.argv[1][:-4]))
    
