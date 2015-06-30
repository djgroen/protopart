# test_io.py

import io 
from domain import *
from nose.tools import assert_true, assert_false  
from nose.tools import assert_raises
import filecmp


def test_ReadDomainGMY():
    d = io.ReadDomain("unittest_data/four_cube.gmy")

#    print d.vertices                #List of all vertices
#    print d.numGraphs             #Number of subgraphs, i.e. number of processors
#    print d.numSites                #Number of sites (Number of vertices per subgraph)
#    print d.neighCount            #Neighbor Count (Number of subgraphs linked to this subgraph)
#    print d.commSurf                #Communication Surface (Number of edges crossing subgraph boundary) 
#    print d.numEdges                #Number Of Edges (Number of edges within a subgraph)
#    print d.totalWeight         #Total Weight according to site Type (Gmy)
#    print d.dictIdx                 #Map vertexID to its index in the vertices list

    assert_true(len(d.vertices) == 64)
#    print len(d.dictIdx)
    assert_true(d.numGraphs == 1)
    assert_true(len(d.dictIdx) == 64)
    
def test_ConvertToCSR():
    d = io.ReadDomain("unittest_data/four_cube.gmy")
#    print len(d.numSites)
#    print d.numGraphs
    csr = io.ConvertToCSR(d)
 
#    print csr.vtxdist
#    print csr.xadj
#    print csr.adjncy
#    print csr.wgts
    
    
def test_WriteCSR():
    d = io.ReadDomain("unittest_data/four_cube.gmy")
    csr = io.ConvertToCSR(d)
    io.WriteCSRFile(csr, "four_cube.csr")
    
def test_ReadCSR():
    c = io.ReadCSRFile("four_cube.csr")
    assert_true(len(c.vtxdist)>0)

def test_WriteCBin():
    d = io.ReadDomain("unittest_data/four_cube.gmy")
    io.writeForCBin(d, "unittest_data/four_cube.hgb")

def test_WriteCBinBif():
    import cProfile, pstats, StringIO
    pr = cProfile.Profile()
    pr.enable()
    d = io.ReadDomain("unittest_data/bif-50um.gmy")
    io.writeForCBin(d, "unittest_data/bif-50um.hgb")
    pr.disable()
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    ps.print_stats()
    print s.getvalue()


def test_ReadAndWriteCBin():
    d = io.readForCBin("unittest_data/four_cube.hgb")
    io.writeForCBin(d, "unittest_data/four_cube2.hgb")
    assert_true(filecmp.cmp("unittest_data/four_cube.hgb", "unittest_data/four_cube2.hgb"))
