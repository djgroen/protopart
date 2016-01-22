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
    
    d = io.ReadDomain(sys.argv[1])

    io.writeHGA(d, "%s.hga" % (sys.argv[1][:-4]))
 
