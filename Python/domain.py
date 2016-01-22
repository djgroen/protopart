import os
import sys
import math
import xdrlib
import zlib
import numpy as np
import itertools
import pymetis
import time
import struct

#Site of the Vertex related to siteType index of VertexGmy object 
siteTypeName = ['bulk','wallsite','inlet','outlet','wallinlet','walloutlet']  
siteTypeWeight = [4,8,16,16,16,16,16]
    
#3D Moore Region coordinate displacement
mooreRegionCor =         [ [-1,-1,-1], [-1,-1, 0], [-1,-1, 1],
                     [-1, 0,-1], [-1, 0, 0], [-1, 0, 1],
                     [-1, 1,-1], [-1, 1, 0], [-1, 1, 1],
                     [ 0,-1,-1], [ 0,-1, 0], [ 0,-1, 1],
                     [ 0, 0,-1],                         [ 0, 0, 1],
                     [ 0, 1,-1], [ 0, 1, 0], [ 0, 1, 1],
                     [ 1,-1,-1], [ 1,-1, 0], [ 1,-1, 1],
                     [ 1, 0,-1], [ 1, 0, 0], [ 1, 0, 1],
                     [ 1, 1,-1], [ 1, 1, 0], [ 1, 1, 1] ]
                     

class Vertex:
    """
    Vertex in graph representation
    """
    def __init__(self):
        self.vertexID = -1                #Vertex ID
        self.coreNum = -1                 #Core Number
        self.edges = []                     #List of Edges it is connected to
        self.edgetypes = []
        self.coordinates = []         #Cartesian Coordinates
        self.siteType = 0
    
class VertexGmy(Vertex):
    """
    Inherited from vertex. Used for Gmy file formats
    due to the extra information
    Note that for GMY Vertices biEdges are bi directional edges to 
    other vertices, while edges are unidirectional edges to 
    a wall site, inlet site, or outlet site
    """
    def __init__(self):
        Vertex.__init__(self)
        self.siteType = -1
        self.biEdges = []
        #self.wallNormal = []         #To be added later

class Domain:
    """
    Graph Class, contains the vertices (vertix Objects)
    """
    
    def __init__(self):
        self.vertices = []             #List of all vertices
        self.numGraphs = 0             #Number of subgraphs, i.e. number of processors
        self.numSites = []             #Number of sites (Number of vertices per subgraph)
        self.neighCount = []         #Neighbor Count (Number of subgraphs linked to this subgraph)
        self.commSurf = []             #Communication Surface (Number of edges crossing subgraph boundary) 
        self.numEdges = []             #Number Of Edges (Number of edges within a subgraph)
        self.totalWeight = []        #Total Weight according to site Type (Gmy)
        self.dictIdx = {}                #Map vertexID to its index in the vertices list
            
    def gmy(self,s):
        """
        Parse GMY file format
        """
        
        #add needed members
        self.gmyFile = file(s)
        self.valid = True
        self.hlbNumber = 0x686c6221
        self.geoMagicNum = 0x676d7904
        self.versionNum = 0x4
        self.blockSize = [] #x,y,z
        self.numLatticeSites = 0 #per one side (^3)
        self.numGraphs = 1
        
        #begin reading/unpacking/decompressing
        self.valid = self._preamble()
        if(not self.valid): return
        self.valid = self._blockHeaders()
        if(not self.valid): return
        self._blockData()
        
            
    def _preamble(self):
        """
        Handles the preamble of the GMY file
        """
        pre = xdrlib.Unpacker(self.gmyFile.read(32))
        
        if (pre.unpack_uint() != self.hlbNumber):
            print 'Incorrect File Peamble!'
            return False
        
        if (pre.unpack_uint() != self.geoMagicNum):
            print 'Incorrect File Peamble!'
            return False
        
        if (pre.unpack_uint() != self.versionNum):
            print 'Incorrect File Peamble!'
            return False
            
        #important information in preamble
        self.blockSize.append(pre.unpack_uint())
        self.blockSize.append(pre.unpack_uint())
        self.blockSize.append(pre.unpack_uint())
        self.numLatticeSites = pre.unpack_uint()

        self.numSites = [[self.numLatticeSites]]
        # GMY files have 1 process only, so numSites on proc 0 must be set equal
        # to the total lattice site count.
        
        if (pre.unpack_uint() != 0):
            print 'Incorrect File Peamble!'
            return False
        
        return True
        
        
    def _blockHeaders(self):
        """
        Saves the information located in the header section
        """
        
        #max range of cartesian coordinates
        self.limit = [self.blockSize[i]*self.numLatticeSites for i in [0,1,2]]
        
        #Header = number of blocks * 3 unsinged integers
        HBytes = self.blockSize[0] * self.blockSize[1] * self.blockSize[2] * 3 * 4
        
        headerLoader = xdrlib.Unpacker(self.gmyFile.read(HBytes))
        
        #number of blocks
        self.nBlocks = self.blockSize[0] * self.blockSize[1] * self.blockSize[2]
    
        #important data
        self.blockFluidSiteCounts = np.zeros(self.nBlocks, dtype=np.uint)
        self.blockDataLength = np.zeros(self.nBlocks, dtype=np.uint)
        self.blockUncompressedDataLength = np.zeros(self.nBlocks, dtype=np.uint)
        #self.blockStarts = np.zeros(self.nBlocks, dtype=np.uint)
        #self.blockEnds = np.zeros(self.nBlocks, dtype=np.uint)
        
        #self.blockEnds[-1] = 32 + HBytes

        #block we are currently in
        self.corBlock = [0,0,-1]
        
        
        for i in list(range(self.nBlocks)):
            #fill the important data
            self.blockFluidSiteCounts[i] = headerLoader.unpack_uint()
            self.blockDataLength[i] = headerLoader.unpack_uint()
            self.blockUncompressedDataLength[i] = headerLoader.unpack_uint()
            #self.blockStarts[i] = self.blockEnds[i-1]
            #self.blockEnds[i] = self.blockStarts[i] + self.blockDataLength[i]
        
	
	self.blockOffset = []
	self.blockIndex = []
	pos = 32 + HBytes
	for i, block_size in enumerate(self.blockDataLength):
		self.blockIndex.append(i)
		self.blockOffset.append(pos)
		pos += int(block_size)


        #Total number of vertices (sites)
        k = 0
        for i in self.blockFluidSiteCounts:
            k += i
        
        self.totalSites = k 
        print 'Total sites: ' + str(k)    
        return True
                
    

                
    def _blockData(self):
        """
        Begin extracting data about sites
        """
        self.dictIdx = {}    #map indeces of vertices
        
        for i in xrange(0, self.nBlocks):
            
            """
            if ((i+1)%80000 == 0):
                print 'Resting for 30sec'
                time.sleep(30)
            """ 
            
            self._updCor(self.corBlock,self.blockSize) #update coordinate of block we are in
            
            if i%1000 == 0:
                print str(i) + '/' + str(self.nBlocks)         #prints block we are reading
                        
            if (self.blockFluidSiteCounts[i] == 0):        #no data
                continue 
            else:
                #decompress, unpack, and handel data
                compressed = self.gmyFile.read(self.blockDataLength[i])
                uncompressed = zlib.decompress(compressed)
                bl = xdrlib.Unpacker(uncompressed)
                self._handleFluid(bl)
                
        
    def _handleFluid(self,fl):
        """
        if a site is fluid, get its edges and add this vertex to list of vertices
        """
        corSite = [0,0,-1] #coordinate within block
        
        TotLatticeSiteCount = self.numLatticeSites**3
        for i in xrange(0, TotLatticeSiteCount):

            self._updCor(corSite, [self.numLatticeSites for k in [0,1,2]]) #increment position within block
            mainType = fl.unpack_uint() 
            if (mainType == 0): #solid, so irrelevant
                continue
            
            #liquid site    
            
            #create Vertex and fill information
            mainV = VertexGmy()
            mainV.coreNum = 0
            xcor = self.corBlock[0] * self.numLatticeSites + corSite[0]
            ycor = self.corBlock[1] * self.numLatticeSites + corSite[1]
            zcor = self.corBlock[2] * self.numLatticeSites + corSite[2]
            mainV.coordinates = [xcor,ycor,zcor]
            mainV.vertexID = self._corToIdx(mainV.coordinates)
            #Error checking
            if self.dictIdx.get(mainV.vertexID, -1) != -1:
                print mainV.vertexID, mainV.coordinates
                x = self.dictIdx[mainV.vertexID]
                print self.vertices[x].vertexID, self.vertices[x].coordinates
            
            self.dictIdx[mainV.vertexID] = len(self.vertices)
           
            for i in xrange(0,26):
                fluidType = fl.unpack_uint()

                #calculate coordinates of Moore Neighbourhood vertex
                edgeCor = [xcor + mooreRegionCor[i][0], ycor + mooreRegionCor[i][1], zcor + mooreRegionCor[i][2]]
                
                
                if (fluidType == 0): #either another vertex or OUT OF BOUNDS
                    mainV.biEdges.append(self._corToIdx(edgeCor))
                    continue              	            
                
                #fluid point (not a vertex) 
                elif (fluidType == 1):
                    #v = VertexWall()
                    #v.distToBoundaryAsFrac
                    x = fl.unpack_float()
                    
                elif (fluidType == 2):
                    #v = VertexInlet()
                    #v.indexToInletArray
                    x = fl.unpack_uint()
                    #v.distToBoundaryAsFrac 
                    x = fl.unpack_float()
                    
                elif (fluidType == 3):
                    #v = VertexOutlet()
                    #v.indexToOutletArray 
                    x = fl.unpack_uint()
                    #v.distToBoundaryAsFrac 
                    x = fl.unpack_float()
                
                
                #append edge with type (Moore, point to a point not a vertex)
                mainV.edges.append(self._corToIdx(edgeCor))
                mainV.edgetypes.append(fluidType)

            #Wall Normal        
            isWallNormal = fl.unpack_uint()
            if (isWallNormal == 1):
                #mainV.wallNormal
                x = [fl.unpack_float(),fl.unpack_float(),fl.unpack_float()] 
            
            
            #Determin Site Type
            if len(mainV.edgetypes) == 0:
                mainV.siteType = 0
            else:
                siteTypeList = list(set(mainV.edgetypes))
                if (siteTypeList == [1]):
                    mainV.siteType = 1
                elif (siteTypeList == [2]):
                    mainV.siteType = 2
                elif (siteTypeList == [3]):
                    mainV.siteType = 3
                elif (siteTypeList == [1,2] or siteTypeList == [2,1]):
                    mainV.siteType = 4
                elif (siteTypeList == [1,3] or siteTypeList == [3,1]) :
                    mainV.siteType = 5
                else:
                    self.siteType = 6 
            
            #add Vertex to Vessel     
            self.vertices.append(mainV)
            
            
    def _updCor(self,cor,L):
        """
        Updates List of 3D coordinates (cor) according to the 
        limitations set by List L
        """
        cor[2] += 1
        if (cor[2] == L[2]):
            cor[1] += 1
            cor[2] = 0
        if (cor[1] == L[1]):
            cor[0] += 1
            cor[1] = 0
        return cor
        
    def _corToIdx(self,L):
        """
        Maps 3D coordinates to a unique index
        """
        return (L[0]*self.limit[1] + L[1])*self.limit[2] + L[2]
        
    def addpoint(self,s):                                             #same
        """
        Format: vertexID, coreNum, coordinates(3), weight, [Edges]
        """
        L = s[0:s.find('[')-1].split()
        
        if (len(L) == 0): 
            return
        
        p = Vertex()
        temp = list(set(s[s.find('[')+1:s.find(']')].split()))
        
        for i in temp: 
            p.edges.append(int(i))
        p.vertexID = int(L[0])
        self.dictIdx[p.vertexID] = len(self.vertices)
        x = p.edges.count(p.vertexID)
        
        if (x > 0):
            p.edges.remove(p.vertexID)
        
        p.coreNum = int(L[1])
        p.coordinates = [float(L[i]) for i in [2,3,4]] 
        
        self.numGraphs = max(self.numGraphs,p.coreNum+1)
        self.vertices.append(p)
        
        
    def changeCoreNum(self,vid,corenum):
        """
        Change the Core Number of a certain vertex by ID
        """
        self.vertices[vid].coreNum = corenum
    
     
    def setNumGraphs(self,x):
        """
        Set the number of subgraphs/partitions
        """
        self.numGraphs = x
        
    
    def disp(self):
        """
        Print Vertices to console with basic information
        """
        for i in self.vertices:
            print 'vertex id: ' + str(i.vertexID)
            print 'core number: ' + str(i.coreNum)
            if (isinstance(i,VertexGmy)):
                print 'edges: ' + ' '.join(str(j) for j in i.biEdges)
            else:
                print 'edges: ' + ' '.join(str(j) for j in i.edges) 
            print 'coordinates: ' + str(i.coordinates)
            print '-------------------------------------------'
        
    def writeo(self,s):
        """
        Print results to file, given as a string containing the path
        """
        wo = file(s,'w')
        for i in self.vertices:
            x = str(i.vertexID) + ' ' + str(i.coreNum) + ' '
            if (not i.coordinates == None):
                x = x + ' '.join(str(j) for j in i.coordinates)
            else: x = x + 'NoCor'
            
            if (isinstance(i,VertexGmy)):
                x += ' SiteType: ' + siteTypeName[i.siteType] + ' ' 
            x = x + ' MooreEdges: ['
            x = x + ' '.join(str(j) for j in i.edges)
            x = x + ']\n'
            if (isinstance(i,VertexGmy)):
                x = x + ' BiEdges: ['
                x = x + ' '.join(str(j) for j in i.biEdges)
                x = x + ']\n' 
            
            wo.write(x)
    
    def updateProperties(self):
        
        for v in self.vertices:
            if v.coreNum > self.numGraphs:
                 self.numGraphs = v.coreNum

        self.numGraphs += 1

        self.numSites = np.zeros(self.numGraphs, dtype=np.uint)
        self.neighCount = []                
        self.commSurf = np.zeros(self.numGraphs, dtype=np.uint)
        self.numEdges = np.zeros(self.numGraphs, dtype=np.uint)
 
    def _getInfo(self):
        """
        Calculates relevant parameters for subgraph division
        """
        self.updateProperties()        
 
        w = {i:[] for i in xrange(0,self.numGraphs)}

        #GMY
        #if (isinstance(self.vertices[0],VertexGmy)):    
        self.totalWeight = [0 for i in list(xrange(0,self.numGraphs))]
        for i in self.vertices:
            if i.coreNum<0 or i.coreNum>self.numGraphs:
                print "Error: coreNum read in is:",i.coreNum
                sys.exit()

            self.numSites[i.coreNum] += 1
            self.totalWeight[i.coreNum] += siteTypeWeight[i.siteType]
            
            edges = i.edges
            if (isinstance(self.vertices[0],VertexGmy)):    
                edges = i.biEdges
                for j in edges:
                    x = self.dictIdx.get(j,-1)
                    if (x == -1): #Out of Bounds
                        continue
                    p = self.vertices[x] 
            
                    if ( (p.coreNum != i.coreNum) 
                           or (j > i.vertexID 
                           and p.coreNum == i.coreNum)):    
                        self.numEdges[i.coreNum] += 1
                    
                    if (p.coreNum != i.coreNum):
                        w[i.coreNum].append(p.coreNum)
                        self.commSurf[i.coreNum] += 1
            else:
                for j in edges:
                    if (j == -1): #Out of Bounds
                        continue
                    p = self.vertices[j]

                    if(p.coreNum == i.coreNum):
                        if j > i.vertexID:
                            self.numEdges[i.coreNum] += 1
                    else:
                        self.numEdges[i.coreNum] += 1
                        w[i.coreNum].append(p.coreNum)
                        self.commSurf[i.coreNum] += 1
   
                    
        for i in w:
            self.neighCount.append(len(set(w[i])))
        
        #prints results of num partitions is small enough
        if (self.numGraphs < 17):
            print 'numSites:'
            print self.numSites
            print 'neighCount:'
            print self.neighCount
            print 'commSurf:'
            print self.commSurf
            print 'numEdges:'
            print self.numEdges
            if (isinstance(self.vertices[0],VertexGmy)):
                print 'totalWeight:'
                print self.totalWeight
    
    def getStats(self):
        self.writeStats(sys.stdout)

    def writeVals(self, f, Avg, maxx, minn):
        f.write( '        Average: ' + str(Avg) + '\n')
        f.write( '        Max        : ' + str(maxx)+ '\n')
        f.write( '        Min        : ' + str(minn)+ '\n')
        if (Avg != 0):
            f.write( '        Max/Avg: ' + str(float(maxx)/Avg) + '\n')

    def writeBlock(self, a):
        summ = 0
        maxx = -1
        minn = a[0]
        for i in a:
            summ += i
            maxx = max(maxx,i)
            minn = min(minn,i)
        Avg = float(summ)/len(a)
        self.writeVals(sys.stdout, Avg, maxx, minn)

    def writeStats(self,f):
        """
        Gives a statistical analysis, writes on file, provided a file
        """
        self._getInfo()
        
        f.write( 'Number of Partitions: ' + str(self.numGraphs) + '\n')
        f.write( '--------------------------------------------\n')
        f.write( ' Number of Sites\n')
        self.writeBlock(self.numSites)
        
        f.write( ' Neighbour Count\n')
        self.writeBlock(self.neighCount)        
 
        f.write( ' Communication Surface\n')
        self.writeBlock(self.commSurf)       
 
        f.write( ' Number of Edges\n')
        self.writeBlock(self.numEdges)       
        
        if (isinstance(self.vertices[0],VertexGmy)):
            f.write( ' Total Weight\n')
            self.writeBlock(self.totalWeight)     
            
        f.write( '--------------------------------------------\n')
