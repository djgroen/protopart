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
    example: python gmy2bin.py <gmy-filename>
    """
    if len(sys.argv) == 1:
        exit(0)
   
    mode = "new"
    if len(sys.argv) > 2:
        mode = sys.argv[2]
 
    if mode == "old":
        d = io.ReadDomain(sys.argv[1])

        d.disp()

        io.writeForCBin(d, "%s.hgb" % (sys.argv[1][:-4]))
    else:
        io.GMY2HGB(sys.argv[1], "%s.hgb" % (sys.argv[1][:-4]), 32)

     
