import os
import sys
import math
import xdrlib
import zlib
import numpy
import itertools
import pymetis
import time
import struct
from domain import *
import io

def splitGraph(d, np, mode="simple"):
    print "Splitting the graph using mode = ", mode
    if mode=="simple":
        splitGraphSimple(d,np)
    elif mode=="recbis":
        splitGraphRecBis(d,np)
    elif mode =="recmul":
        splitGraphRecMul(d,np)
    elif mode =="recbisweighted":
        splitGraphRecMulWeighted(d,np)
    elif mode =="recmulweighted":
        splitGraphRecMulWeighted(d,np)
    elif mode =="metis":
        splitGraphMetis(d,np)
    elif mode =="metisweighted":
        splitGraphMetisWeighted(d,np)
    elif mode =="recmetis":
        splitGraphRecMetis(d,np)
    elif mode =="recmetisweighted":
        splitGraphRecMetisWeighted(d,np)

def splitGraphSimple(d,np):
    """
    Split to subgraphs by dividing number of sites and number of 
    edges evenly among subgraphs (coreNum)
    """
    d.numGraphs = np
    count = 0
    for i in d.vertices:
        if (isinstance(i,VertexGmy)):
            count += 1 + len(i.biEdges)
        else:
            count = count + 1 + len(i.edges)
    count = count/np
    
    proc = 0
    j = 0
        
    for i in d.vertices:
        i.coreNum = proc;
        if (isinstance(i,VertexGmy)):
            j += 1 + len(i.biEdges)
        else:
            j = j + 1 + len(i.edges)
        if (j >= count):
            j = 0
            proc = proc+1
    
    
    
def splitGraphRecBis(d, np):
    """
    Split to subgraphs using the Recursive Bisection Method. 
    Cartesian Coordinates Required. Must be divided into powers of 2
    """
    if (np & (np-1) != 0): 
        print 'Error, not power of 2!'
        return
        
    d.numGraphs = np
    for i in d.vertices: #initialize all core numbers to zero
        i.coreNum = 0
    np = int(math.log(np,2))
    W = list(range(0,len(d.vertices)))
    _splitGraphRecBisPriv(d, W, np-1)                 
        
def _splitGraphRecBisPriv(d, W, call):
    if (call == -1): return
    L = []
    for i in W:
        L.append( (d.vertices[i].coordinates[call%3],i) )
    L.sort()
    
    for i in L[len(L)/2:]:
        d.vertices[i[1]].coreNum += 2**call
            
    _splitGraphRecBisPriv(d, [i[1] for i in L[:len(L)/2]], call-1)
    _splitGraphRecBisPriv(d, [i[1] for i in L[len(L)/2:]], call-1)

def _primeFac(n):
    """
    Returns a list of all prime factors of n
    """
    step = lambda x: 1 + x*4 - (x//2)*2
    maxq = long(math.floor(math.sqrt(n)))
    d = 1
    q = n % 2 == 0 and 2 or 3 
    while q <= maxq and n % q != 0:
        q = step(d)
        d += 1
    res = []
    if q <= maxq:
        res.extend(_primeFac(n//q))
        res.extend(_primeFac(q)) 
    else: res=[n]
    return res

def splitGraphRecMul(d, x):
    """
    Split to subgraphs using the Recursive Multisection Method. 
    Cartesian Coordinates Required.
    """
    d.numGraphs = x
    for i in d.vertices: #initialize all core numbers to zero
        i.coreNum = 0
    
    W = list(range(0,len(d.vertices)))
    pf = _primeFac(x)
    _splitGraphRecMulPriv(d, W, x, len(pf)-1, pf)                    


def _splitGraphRecMulPriv(d, W, call, idx, pf):
    
    if (idx == -1): return
    L = []
    for i in W:
        L.append( (d.vertices[i].coordinates[idx%3],i) )
    L.sort()
        
        
    for i in list(range(pf[idx])):
        Q = L[i*len(L)/pf[idx]:(i+1)*len(L)/pf[idx]]
        for j in Q:
            d.vertices[j[1]].coreNum += i*(call/pf[idx])
        _splitGraphRecMulPriv(d, [q[1] for q in Q],call/pf[idx],idx-1,pf)
        
def splitGraphRecMulWeighted(d, x):
    """
    Split to subgraphs using the Recursive Multisection Method. 
    Cartesian Coordinates Required. Uses Site Weights
    """
    d.numGraphs = x
    for i in d.vertices: #initialize all core numbers to zero
        i.coreNum = 0
    
    W = list(range(0,len(d.vertices)))
    pf = _primeFac(x)
    _splitGraphRecMulPrivWeighted(d, W,x,len(pf)-1,pf)                    


def _splitGraphRecMulPrivWeighted(d, W, call, idx, pf):
    
    if (idx == -1): return
    L = []
    for i in W:
        L.append( (d.vertices[i].coordinates[idx%3],i,siteTypeWeight[d.vertices[i].siteType]) )
    L.sort()
    
    # Weighting
    weightSum = 0
    for i in L:
        weightSum += i[2]
    target = weightSum/pf[idx]
    
    indexList = [0]
    weightSum = 0
    for i in list(range(len(L))):
        weightSum += L[i][2]
        if (weightSum >= target):
            indexList.append(i+1)
            weightSum = 0
    
    if (indexList[-1] != len(L)):
        indexList.append(len(L))        
    
    #
    
    
    for i in list(range(len(indexList)-1)):
        Q = L[indexList[i]:indexList[i+1]]
        for j in Q:
            d.vertices[j[1]].coreNum += i*(call/pf[idx])
        _splitGraphRecMulPrivWeighted(d, [q[1] for q in Q],call/pf[idx],idx-1,pf)



def splitGraphRecBisWeighted(d, x):
    """
    Split to subgraphs using the Recursive Bisection Method. 
    Cartesian Coordinates Required. Must be divided into powers of 2
    Uses the provided Weights
    """
    if (x & (x-1) != 0): 
        print 'Error, not power of 2!'
        return
        
    d.numGraphs = x
    for i in d.vertices: #initialize all core numbers to zero
        i.coreNum = 0
    x = int(math.log(x,2))
    W = list(range(0,len(d.vertices)))
    _splitGraphRecBisPrivWeighted(W,x-1)                 
        

def _splitGraphRecBisPrivWeighted(d, W, call):
    if (call == -1): return
    L = []
    for i in W:
        L.append( (d.vertices[i].coordinates[call%3],i,siteTypeWeight[d.vertices[i].siteType]) )
    L.sort()
    
    # Weighting
    weightSum = 0
    for i in L:
        weightSum += i[2]
    target = weightSum/2
    
    weightSum = 0
    idx = 0
    for i in list(range(len(L))):
        weightSum +=    L[i][2]
        if (weightSum >= target):
            idx = i+1
            break
    #
    for i in L[idx:]:
        d.vertices[i[1]].coreNum += 2**call
            
    _splitGraphRecBisPrivWeighted(d, [i[1] for i in L[:idx]], call-1)
    _splitGraphRecBisPrivWeighted(d, [i[1] for i in L[idx:]], call-1)
    

def splitGraphMetis(d,x):
    """
    Splits the graph using Metis library function, iterative, no weights
    """
    L = {}
    if (isinstance(d.vertices[0],VertexGmy)):
        for i in d.vertices:
            q = []
            for j in i.biEdges:
                m = d.dictIdx.get(j,-1)
                if (m != -1):
                    q.append(m)
            L[d.dictIdx[i.vertexID]]= q
    else:
        for i in d.vertices:
            q = []
            for j in i.edges:
                m = d.dictIdx.get(j,-1)
                if (m != -1):
                    q.append(m)
            L[d.dictIdx[i.vertexID]]= q
    
    print 'Done Making Graph!'
    
    (wut,edges) = pymetis.part_graph(x,L)
    
    #Determine the number of splits, since Metis
    #sometimes returns less than the number of 
    #partitions requested by user
    numsplits = []
    for i in edges:
        numsplits.append(i)
    numsplits = len(list(set(numsplits)))
    
    if numsplits != x:
        print 'Numsplits = ' + str(numsplits) + ', not ' + str(x) + '!'
        
    j = 0
    for i in d.vertices: 
        i.coreNum = edges[j]
        j += 1
    
    d.numGraphs = numsplits
    

def splitGraphRecMetis(d,x):
    """
    Splits the graph using Metis library function, recursive, no weights
    """
    L = {}
    if (isinstance(d.vertices[0],VertexGmy)):
        for i in d.vertices:
            q = []
            for j in i.biEdges:
                m = d.dictIdx.get(j,-1)
                if (m != -1):
                    q.append(m)
            L[d.dictIdx[i.vertexID]]= q
    else:
        for i in d.vertices:
            q = []
            for j in i.edges:
                m = d.dictIdx.get(j,-1)
                if (m != -1):
                    q.append(m)
            L[d.dictIdx[i.vertexID]]= q
    
    print 'Done Making Graph!'
    
    (wut,edges) = pymetis.part_graph(x,L,None,None,None,None,True)
    
    #Determine the number of splits, since Metis
    #sometimes returns less than the number of 
    #partitions requested by user
    numsplits = []
    for i in edges:
        numsplits.append(i)
    numsplits = len(list(set(numsplits)))
    
    if numsplits != x:
        print 'Numsplits = ' + str(numsplits) + ', not ' + str(x) + '!'
        
    j = 0
    for i in d.vertices: 
        i.coreNum = edges[j]
        j += 1
    
    d.numGraphs = numsplits


def splitGraphMetisWeighted(d,x):
    """
    Splits the graph using Metis library function, iterative, with weights
    """
    L = {}
    wght = []
    for i in d.vertices:
        q = []
        wght.append(siteTypeWeight[i.siteType])
        for j in i.biEdges:
            m = d.dictIdx.get(j,-1)
            if (m != -1):
                q.append(m)
        L[d.dictIdx[i.vertexID]]=q
    
    print 'Done Making Graph!'
    
    (wut,edges) = pymetis.part_graph(x,L,None,None,wght,None)
    
    #Determine the number of splits, since Metis
    #sometimes returns less than the number of 
    #partitions requested by user
    numsplits = []
    for i in edges:
        numsplits.append(i)
    numsplits = len(list(set(numsplits)))
    
    if numsplits != x:
        print 'Numsplits = ' + str(numsplits) + ', not ' + str(x) + '!'
        
    j = 0
    for i in d.vertices: 
        i.coreNum = edges[j]
        j += 1
    
    d.numGraphs = numsplits
    
    
def splitGraphRecMetisWeighted(d,x):
    """
    Splits the graph using Metis library function, recursive, no with weights
    """
    L = {}
    wght = []
    for i in d.vertices:
        q = []
        wght.append(siteTypeWeight[i.siteType])
        for j in i.biEdges:
            m = d.dictIdx.get(j,-1)
            if (m != -1):
                q.append(m)
        L[d.dictIdx[i.vertexID]]= q
    
    
    print 'Done Making Graph!'
    
    (wut,edges) = pymetis.part_graph(x,L,None,None,wght,None,True)
    
    numsplits = []
    for i in edges:
        numsplits.append(i)
    numsplits = len(list(set(numsplits)))
    
    if numsplits != x:
        print 'Numsplits = ' + str(numsplits) + ', not ' + str(x) + '!'
        
    j = 0
    for i in d.vertices: 
        i.coreNum = edges[j]
        j += 1
    
    d.numGraphs = numsplits

