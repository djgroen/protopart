# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 12:05:02 2014

@author: derek
"""
import io
import decomposition as dec
import sys
import domain


def mainall(W,x):

    v = io.ReadDomain(W['input'])
    
    wout = open(W['output'],'w')
    
    if (isinstance(v.vertices[0],domain.VertexGmy)):
        wout.write("Block Size: " + str(v.blockSize)+ '\n')
        wout.write("Num Lattice Sites per site per block: " + str(v.numLatticeSites) + '\n')
        wout.write("TotalSites: " + str(v.totalSites) + '\n')
        wout.write( '--------------------------------------------\n')
        wout.write( '--------------------------------------------\n')
        wout.write( '--------------------------------------------\n')

    modes = ["simple","recbis","recmul","recbisweighted","recmulweighted"]    
    for m in modes:
        # exclude "metis", "metisweighted","recmetis", "recmetisweighted" for now.
        dec.splitGraph(v, x, mode=m)
        v.writeStats(wout)
                    
    #if (isinstance(v.vertices[0],VertexGmy)):
    #This statements detects whether the source is a gmy file.
        
    if (W['outputc'] != ''):
        io.writeHGA(v, W['outputc'])
    
    if (W['outputb'] != ''):
        io.writeForCBin(v, W['outputb'])
    
    
def mainone(W,x):
    
    v = io.ReadDomain(W['input'])
        
    func = ''

    print 'Using mode ', W['decomp']
    dec.splitGraph(x, mode=W['decomp'])
    
    if (W['output'] != ''): 
        wout = open(W['output'],'w')
        
        if (isinstance(v.vertices[0],domain.VertexGmy)):
            wout.write("Block Size: " + str(v.blockSize)+ '\n')
            wout.write("Num Lattice Sites per site per block: " + str(v.numLatticeSites) + '\n')
            wout.write("TotalSites: " + str(v.totalSites) + '\n')
            wout.write( '--------------------------------------------\n')
        wout.write(func+'\n')
        v.writeStats(wout)

    if (W['outputc'] != ''):
        io.writeHGA(v, W['outputc'])
    
    if (W['outputb'] != ''):
        io.writeForCBin(v, W['outputb'])
        
        
if __name__ == '__main__':
    """
    -i input file
    -o output file for decomposition quality
    -c output file for c++ as txt
    -b output file for c++ as bin
        
    --all : run all decomposition methods
    if only one:
    --mode <name of mode>
    
    integer: num partitions
    
    example: python decomposition.py -i input.gmy -o output.txt -b outputc.bin --all 20
    """
    cmdargs = list(sys.argv)
    if len(cmdargs) == 1:
        exit(0)
    
    D = {'input': '', 'output':'', 'outputc':'', 'outputb':'','decomp':''}
    
    numPartitions = int(cmdargs[-1])
    
    for i in list(range(len(cmdargs)-1)):
        if (cmdargs[i][0] != '-'):
            continue
        if (cmdargs[i] == '-i'):
            D['input'] = cmdargs[i+1]
        if (cmdargs[i] == '-o'):
            D['output'] = cmdargs[i+1]
        if (cmdargs[i] == '-c'):
            D['outputc'] = cmdargs[i+1]
        if (cmdargs[i] == '-b'):
            D['outputb'] = cmdargs[i+1]
        if (cmdargs[i] == '--all'):
            D['decomp'] = 'all'
        if (cmdargs[i] == '--mode'):
            D['decomp'] = cmdargs[i+1]
            
    if (D['decomp'][0:4] == 'all'):
        mainall(D,numPartitions)
    else:
        mainone(D,numPartitions)
