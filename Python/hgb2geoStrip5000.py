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

    d = io.readForCBin(sys.argv[1])

    print "%s.gmyhgb" % (sys.argv[1][:-5])
    d2 = io.ReadDomain("%s.gmy" % (sys.argv[1][:-5]))
	
    d.numLatticeSites = d2.numLatticeSites
    d.blockSize = d2.blockSize
    d.blockFluidSiteCounts = d2.blockFluidSiteCounts
    d.blockOffset = d2.blockOffset
    d.blockDataLength = d2.blockDataLength
    d.blockUncompressedDataLength = d2.blockUncompressedDataLength

    #d.disp()

    io.writeForCBinX35_987(d, "%s.gmyhgb" % (sys.argv[1][:-5]))
    io.writeForCBinPlus4000(d, "%s.gmyrnk" % (sys.argv[1][:-5]))
    
