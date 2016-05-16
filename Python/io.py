import domain
import csr
import os
import numpy as np
import struct  
import xdrlib
import math    
import zlib
import xdrlib

def extractPropertiesFromFluid(d,fl,mode="neighbour"):
    """
    if a site is fluid, get its edges and add this vertex to list of vertices
    Mode = neighbour: add and return only neighbour info
    Mode = neighbour: add and return only d info
    """
    
    headers = [[],[],[]] #create a list of headers for "neighbour" mode. This will contain [[vertexID],[dictIdx],[biEdges]]
    vertices = [] #create a list of vertices for "self" mode.
    
    corSite = [0,0,-1] #coordinate within block
    
    elements_passed = 0

    TotLatticeSiteCount = d.numLatticeSites**3
    for i in xrange(0, TotLatticeSiteCount):

        d._updCor(corSite, [d.numLatticeSites for k in [0,1,2]]) #increment position within block
        mainType = fl.unpack_uint() 
        if (mainType == 0): #solid, so irrelevant
            continue
        
        #liquid site    
        
        #create Vertex and fill information
        mainV = domain.VertexGmy()
        
        xcor = d.corBlock[0] * d.numLatticeSites + corSite[0]
        ycor = d.corBlock[1] * d.numLatticeSites + corSite[1]
        zcor = d.corBlock[2] * d.numLatticeSites + corSite[2]
        mainV.coordinates = [xcor,ycor,zcor]
        mainV.vertexID = d._corToIdx(mainV.coordinates)
        mainV.coreNum = mainV.vertexID / d.numGraphs # numGraphs is the number of cores planned.
        #Error checking
        if d.dictIdx.get(mainV.vertexID, -1) != -1:
            #print mainV.vertexID, mainV.coordinates
            x = d.dictIdx[mainV.vertexID]
            #print d.vertices[x].vertexID, d.vertices[x].coordinates
        
        if mode == "neighbour":
            headers[0].append(mainV.vertexID)
            headers[1].append(elements_passed)
            headers[2].append([])
        
        ### dictIdx STILL GROWS with NumSites.
        d.dictIdx[mainV.vertexID] = elements_passed
        elements_passed += 1
            
        for i in xrange(0,26):
            fluidType = fl.unpack_uint()

            #calculate coordinates of Moore Neighbourhood vertex
            edgeCor = [xcor + domain.mooreRegionCor[i][0], ycor + domain.mooreRegionCor[i][1], zcor + domain.mooreRegionCor[i][2]]
            
            
            if (fluidType == 0): #either another vertex or OUT OF BOUNDS
                mainV.biEdges.append(d._corToIdx(edgeCor))
                if mode == "neighbour":
                    headers[2][-1].append(d._corToIdx(edgeCor))
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
            mainV.edges.append(d._corToIdx(edgeCor))
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
                mainV.siteType = -1 
        
        #add Vertex to Vessel     
        vertices.append(mainV)
        
        
    if mode == "self":
        return vertices
    if mode == "neighbour":
        return headers




def GMY2HGB(gmyfile, hgbfile, corecount, disp=True):

    f = open(hgbfile,'wb') #open output HGB File

    ### READING PASS 1 (header)
    d = domain.Domain()
    
    #Check if file exists
    if (not os.path.isfile(gmyfile)):     
        print 'Error! Input file not found!'
        return
        
    
    #Determine File Format Type
    formatType = 0; 
   
    d.gmy(gmyfile, readData=False) #Read in header information from the GMY, but ignore the data fields.
   

    ### WRITING PART
    # 1. Total site aocunt
    x = struct.pack('i',np.sum(d.blockFluidSiteCounts))
    print "Site count = ", np.sum(d.blockFluidSiteCounts)
    f.write(x)    
    
    
    ### READING PASS 1 (fluid sites)
    d.dictIdx = {}    #map indices of vertices

    neighbour_data = []
    
    for i in xrange(0, d.nBlocks):

        """
        if ((i+1)%80000 == 0):
            print 'Resting for 30sec'
            time.sleep(30)
        """

        d.updCor() #update coordinate of block we are in

        if i%1000 == 0:
            print str(i) + '/' + str(d.nBlocks)         #prints block we are reading

        if (d.blockFluidSiteCounts[i] == 0):        #no data
            continue
        else:
            #decompress, unpack, and handel data
            compressed = d.gmyFile.read(d.blockDataLength[i])
            uncompressed = zlib.decompress(compressed)
            bl = xdrlib.Unpacker(uncompressed)
            neighbour_data = extractPropertiesFromFluid(d,bl,mode="neighbour")

            # WRITING (neighbour information)
            # 2. neighbour counts for each Vertex
            for i in xrange(0,len(neighbour_data[0])):
                # x = struct.pack('i',i.vertexID)
                # print i.vertexID, d.dictIdx[i.vertexID]
                x = struct.pack('i',neighbour_data[1][i])
                f.write(x)
                
                x = struct.pack('i',len(neighbour_data[2][i]))
                f.write(x)
                # print "writing: ", d.dictIdx[i.vertexID], len(i.biEdges)

                
    ### READING PASS 2 (Ignore the headers)
    d = domain.Domain()
   
    d.gmy(gmyfile, readData=False) #Read in header information from the GMY, but ignore the data fields.
    d.setNumGraphs(corecount)   
 
    ### READING PASS 2 (fluid sites)
    vertices = []
    
    for i in xrange(0, d.nBlocks):

        """
        if ((i+1)%80000 == 0):
            print 'Resting for 30sec'
            time.sleep(30)
        """

        d.updCor() #update coordinate of block we are in

        if i%1000 == 0:
            print str(i) + '/' + str(d.nBlocks)         #prints block we are reading

        if (d.blockFluidSiteCounts[i] == 0):        #no data
            continue
        else:
            #decompress, unpack, and handel data
            compressed = d.gmyFile.read(d.blockDataLength[i])
            uncompressed = zlib.decompress(compressed)
            bl = xdrlib.Unpacker(uncompressed)
            vertices = extractPropertiesFromFluid(d,bl,mode="self")
 
            # WRITING: Vertex data 
            for i in vertices:
                print 'vertex id: ' + str(d.dictIdx[i.vertexID])
                print 'core number: ' + str(i.coreNum)
                print 'site type: ' + str(i.siteType), i.edgetypes
                print 'edges: ' + ' '.join(str(j) for j in i.biEdges)
                print 'coordinates: ' + str(i.coordinates)
                print '-------------------------------------------'


                # x = struct.pack('i',i.vertexID)
                x = struct.pack('i',d.dictIdx[i.vertexID])
                f.write(x)
                x = struct.pack('i',int( (float(d.dictIdx[i.vertexID]) / float(len(d.dictIdx)) ) * d.numGraphs))
                f.write(x)
                x = struct.pack('i',i.coordinates[0])
                f.write(x)
                x = struct.pack('i',i.coordinates[1])
                f.write(x)
                x = struct.pack('i',i.coordinates[2])
                f.write(x)
                x = struct.pack('i',i.siteType)
                f.write(x)
        
                for j in i.biEdges:
                    x = d.dictIdx.get(j,-1)
                    if (x == -1): continue
                    x = struct.pack('i',d.dictIdx[j])
                    f.write(x)
        
                b = []
                for j in i.biEdges:
                    b.append(d.dictIdx.get(j,-1))    
                # print "writing: ", d.dictIdx[i.vertexID], i.coreNum, i.coordinates, i.siteType, b



def ReadDomain(rfile):
    """
    Takes as input string containing input file path, determines
    the format, and calls the appropriate parse function
    """
        
    v = domain.Domain()
    
    #Check if file exists
    if (not os.path.isfile(rfile)):     
        print 'Error! Input file not found!'
        return
        
    
    #Determine File Format Type
    formatType = -1;
    
    if (rfile[-4:] == '.gmy'):
        formatType = 0
    else:     
        xfile = file(rfile,'r') 
        s = xfile.readline()
    
        formatType = len(s[0:s.find('[')-1].split())
    
    #Format Type 1: vertexID, coreNum, coordinates(3), weights, [Edges]
    #Format Type 0: .gmy file
    
    
    #call function according to format Type
    if (formatType == 0):
        v.gmy(rfile)
    elif (formatType == 1):
        while (not s == ''):
            v.addpoint(s)
            s = xfile.readline()
    
    return v


def writeHGA(d, path, weights_only=False):
    """
    Writes a file for c++ to parse in ASCII format, provided a path
    Core ID - Coordinates - SiteTypeNum - Edges - \n
    """
    print "Writing HGA file at: ", path
    f = open(path,'w')

    s = str(len(d.vertices)) + "\n"

    if weights_only:
        s = ""
        for i in d.vertices:
            if i.siteType == 2 or i.siteType == 4:
                s += str(i.coordinates[0]) + ' '
                s += str(i.coordinates[1]) + ' '
                s += str(i.coordinates[2]) + ' 0.0\n'
        f.write(s)
        return


    for i in d.vertices:
        s += str(d.dictIdx[i.vertexID]) + ' ' + str(len(i.biEdges)) + "\n"

    for i in d.vertices:
        s += str(d.dictIdx[i.vertexID]) + ' ' + str(i.coreNum) + ' ' + str(i.coordinates[0]) + ' '
        s += str(i.coordinates[1]) + ' ' + str(i.coordinates[2]) + ' '
        s += str(i.siteType) + ' '
        
        for j in i.biEdges:
            x = d.dictIdx.get(j,-1)
            if (x == -1): continue
            s += str(d.dictIdx[j]) + ' '
        s += '\n'

    f.write(s)

def writeVelWeights(d, path):
    """
    Writes a file for c++ to parse in ASCII format, provided a path
    Coordinates - Weight (set to 0.0) \n
    """
    print "Writing HGA file at: ", path
    f = open(path,'w')

    for i in d.vertices:
        s = str(i.coordinates[0]) + ' ' + str(i.coordinates[1]) + ' ' + str(i.coordinates[2]) + ' 0.0'
        s += str(i.siteType) + ' '
        s += '\n'

        f.write(s)

def writeHGAX35_987(d, path):
    
    # X35_987 is the XDR version of the output

    """
    Writes a file for c++ to parse in binary format, provided a path
    Core ID - Coordinates - SiteTypeNum - Edges - -1
    """
    f = open(path,'w')
    entries = []

    # 1. Total site aocunt
    # struct.pack('i',totalSites)
    print "Site count = ", len(d.vertices)

    blockcoord = [0.0,0.0,0.0]

    # 3. Vertex data 
    for i in d.vertices:

# Read halo- blocks on domain edges
      if (i.coordinates[0]%d.numLatticeSites == 0): #Lower x-edge
        xRange = [-1,0]
      elif (i.coordinates[0]%d.numLatticeSites == 7): #Upper x-edge
        xRange = [0,1]
      else:
        xRange = [0]

      if (i.coordinates[1]%d.numLatticeSites == 0): #Lower x-edge
        yRange = [-1,0]
      elif (i.coordinates[1]%d.numLatticeSites == 7): #Upper x-edge
        yRange = [0,1]
      else:
        yRange = [0]

      if (i.coordinates[2]%d.numLatticeSites == 0): #Lower x-edge
        zRange = [-1,0]
      elif (i.coordinates[2]%d.numLatticeSites == 7): #Upper x-edge
        zRange = [0,1]
      else:
        zRange = [0]


      for dx in xRange:
        for dy in yRange:
          for dz in zRange:
            blockcoord[0] = math.floor((i.coordinates[0]+dx)/d.numLatticeSites)
            blockcoord[1] = math.floor((i.coordinates[1]+dy)/d.numLatticeSites)
            blockcoord[2] = math.floor((i.coordinates[2]+dz)/d.numLatticeSites)
            blockid=int((blockcoord[0]*d.blockSize[1]+blockcoord[1])*d.blockSize[2]+blockcoord[2])

            key = (i.coreNum,blockid)
            if key not in entries:
                entries.append(key)

    for i in d.vertices:
      blockcoord[0] = math.floor((i.coordinates[0])/d.numLatticeSites)
      blockcoord[1] = math.floor((i.coordinates[1])/d.numLatticeSites)
      blockcoord[2] = math.floor((i.coordinates[2])/d.numLatticeSites)
      blockID=int((blockcoord[0]*d.blockSize[1]+blockcoord[1])*d.blockSize[2]+blockcoord[2])

      siteID=int((i.coordinates[0]%d.numLatticeSites)*d.numLatticeSites+(i.coordinates[1]%d.numLatticeSites))*d.numLatticeSites+(i.coordinates[2]%d.numLatticeSites)
      cpus = []

      for key in entries:
        core = key[0]
        block = key[1]
        if (block==blockID):
          cpus.append(core)

      print d.dictIdx[i.vertexID],i.coreNum,i.coordinates[0],i.coordinates[1],i.coordinates[2],i.siteType,len(i.biEdges)
 
      for j in i.biEdges:
          x = d.dictIdx.get(j,-1)
          if (x == -1): continue
          print d.dictIdx[j]
  
      print blockID,d.blockFluidSiteCounts[blockID],siteID,len(cpus)

      for j in cpus:
        print j

      s=str(i.coordinates[0]) + ' ' + str(i.coordinates[1]) + ' ' + str(i.coordinates[2]) + ' ' + str(i.siteType) + ' ' + str(i.coreNum) + ' ' + str(blockID) + ' '

      for j in cpus:
        s+= str(j) + ' '

      s+="\n"

      f.write(s)


def writeHGAplus4000(d, path):
    """
    Writes a file for c++ to parse in ASCII format, provided a path
    Core ID - Coordinates - SiteTypeNum - Edges - \n
    """
    entries = []
    print "Writing HGA file at: ", path
    f = open(path,'w')
    blockcoord = [0.0,0.0,0.0]
    for i in d.vertices:
    # Include neighbouring blocks when on domain edges
      if (i.coordinates[0]%d.numLatticeSites == 0): #Lower x-edge
        xRange = [-1,0]
      elif (i.coordinates[0]%d.numLatticeSites == 7): #Upper x-edge
        xRange = [0,1]
      else:
        xRange = [0]

      if (i.coordinates[1]%d.numLatticeSites == 0): #Lower x-edge
        yRange = [-1,0]
      elif (i.coordinates[1]%d.numLatticeSites == 7): #Upper x-edge
        yRange = [0,1]
      else:
        yRange = [0]

      if (i.coordinates[2]%d.numLatticeSites == 0): #Lower x-edge
        zRange = [-1,0]
      elif (i.coordinates[2]%d.numLatticeSites == 7): #Upper x-edge
        zRange = [0,1]
      else:
        zRange = [0]
      
      for dx in xRange:
        for dy in yRange:
          for dz in zRange:
            blockcoord[0] = math.floor((i.coordinates[0]+dx)/d.numLatticeSites)
	    blockcoord[1] = math.floor((i.coordinates[1]+dy)/d.numLatticeSites)
            blockcoord[2] = math.floor((i.coordinates[2]+dz)/d.numLatticeSites)
	    # Recalculate the blockcoords
	    blockid=int((blockcoord[0]*d.blockSize[1]+blockcoord[1])*d.blockSize[2]+blockcoord[2])

	    key = (i.coreNum,blockid)
            if key not in entries:
		entries.append(key)
	       
                s = str(i.coreNum) + ' ' + str(blockid) + ' ' + str(d.blockOffset[blockid]) + ' ' + str(d.blockDataLength[blockid]) + ' ' + str(d.blockUncompressedDataLength[blockid]) + ' ' + str(d.blockFluidSiteCounts[blockid])

                s += '\n'
                f.write(s)



def readHGA(path):
    """
    Writes a file for c++ to parse in ASCII format, provided a path
    Core ID - Coordinates - SiteTypeNum - Edges - \n
    """
    print "Reading HGA file at: ", path

    d = domain.Domain()

    f = open(path,'r')
 
    #lines = f.readlines()

    lineCount = 0
    for line in f:
       s = (line.strip()).split()
       v = domain.Vertex()
       v.coreNum = int(s[0])
       v.coordinates = np.array([int(s[1]), int(s[2]), int(s[3])])
       v.siteType = int(s[4])
       for i in xrange(5, len(s)):
           v.edges.append(int(s[i]))
       lineCount += 1
       d.vertices.append(v)

    return d
"""
    for i in d.vertices:
        s = str(i.coreNum) + ' ' + str(i.coordinates[0]) + ' '
        s += str(i.coordinates[1]) + ' ' + str(i.coordinates[2]) + ' '
        s += str(i.siteType) + ' '

        for j in i.biEdges:
            x = d.dictIdx.get(j,-1)
            if (x == -1): continue
            s += str(d.dictIdx[j]) + ' '
        s += '\n'
        f.write(s)

    d.updateProperties()
    return d
"""
        
def writeMapping(d, path):
    """
    Writes a file in ASCII format, provided a path
    Core ID - Coordinates \n
    """
    print "Writing Mapping file at: ", path
    f = open(path,'w')
    for i in d.vertices:
        s = str(i.coreNum) + ' ' + str(i.coordinates[0]) + ' '
        s += str(i.coordinates[1]) + ' ' + str(i.coordinates[2]) + ' '
        s += '\n'
        f.write(s)

def writeForCBin(d, path):
    """
    Writes a file for c++ to parse in binary format, provided a path
    Core ID - Coordinates - SiteTypeNum - Edges - -1
    """
    f = open(path,'wb')

    # 1. Total site aocunt
    # x = struct.pack('i',totalSites)
    x = struct.pack('i',len(d.vertices))
    print "Site count = ", len(d.vertices)
    f.write(x)    

    # 2. neighbour counts for each Vertex
    for i in d.vertices:
        # x = struct.pack('i',i.vertexID)
        # print i.vertexID, d.dictIdx[i.vertexID]
        x = struct.pack('i',d.dictIdx[i.vertexID])
 
        f.write(x)
        x = struct.pack('i',len(i.biEdges))
        f.write(x)
        # print "writing: ", d.dictIdx[i.vertexID], len(i.biEdges)

 
    # 3. Vertex data 
    for i in d.vertices:
#        x = struct.pack('i',i.vertexID)
        x = struct.pack('i',d.dictIdx[i.vertexID])
        f.write(x)
        x = struct.pack('i',i.coreNum)
        f.write(x)
        x = struct.pack('i',i.coordinates[0])
        f.write(x)
        x = struct.pack('i',i.coordinates[1])
        f.write(x)
        x = struct.pack('i',i.coordinates[2])
        f.write(x)
        x = struct.pack('i',i.siteType)
        f.write(x)
        
        for j in i.biEdges:
            x = d.dictIdx.get(j,-1)
            if (x == -1): continue
            x = struct.pack('i',d.dictIdx[j])
            f.write(x)
        
        b = []
        for j in i.biEdges:
            b.append(d.dictIdx.get(j,-1))    
        # print "writing: ", d.dictIdx[i.vertexID], i.coreNum, i.coordinates, i.siteType, b


def populateBlockList(v, d):
    """
    Creates a mapping list from blocknumber to corenumber.
    """
    entries = []

    #read halo- blocks on domain edges
    if (i.coordinates[0]%d.numLatticeSites == 0): #Lower x-edge
        xRange = [-1,0]
    elif (i.coordinates[0]%d.numLatticeSites == 7): #Upper x-edge
        xRange = [0,1]
    else:
        xRange = [0]

    if (i.coordinates[1]%d.numLatticeSites == 0): #Lower x-edge
        yRange = [-1,0]
    elif (i.coordinates[1]%d.numLatticeSites == 7): #Upper x-edge
        yRange = [0,1]
    else:
        yRange = [0]

    if (i.coordinates[2]%d.numLatticeSites == 0): #Lower x-edge
        zRange = [-1,0]
    elif (i.coordinates[2]%d.numLatticeSites == 7): #Upper x-edge
        zRange = [0,1]
    else:
        zRange = [0]

    blockcoord = [0.0,0.0,0.0]
    for dx in xRange:
      for dy in yRange:
        for dz in zRange:
          blockcoord[0] = math.floor((i.coordinates[0]+dx)/d.numLatticeSites)
          blockcoord[1] = math.floor((i.coordinates[1]+dy)/d.numLatticeSites)
          blockcoord[2] = math.floor((i.coordinates[2]+dz)/d.numLatticeSites)
          blockid=int((blockcoord[0]*d.blockSize[1]+blockcoord[1])*d.blockSize[2]+blockcoord[2])

          key = (i.coreNum,blockid)
          if key not in entries:
            entries.append(key)



def writeForCBinX35_987(d, path):
    
    # X35_987 is the XDR version of the output

    """
    Writes a file for c++ to parse in binary format, provided a path
    Core ID - Coordinates - SiteTypeNum - Edges - -1
    """
    f = open(path,'wb')

    p = xdrlib.Packer()

    # 1. Total site aocunt
    # struct.pack('i',totalSites)
    p.pack_uint(len(d.vertices))
    print "Site count = ", len(d.vertices)

    blockcoord = [0.0,0.0,0.0]

    # 3. Vertex data 
    for i in d.vertices:
      entries = populateBlockList(i, d)

    for i in d.vertices:
      blockcoord[0] = math.floor((i.coordinates[0])/d.numLatticeSites)
      blockcoord[1] = math.floor((i.coordinates[1])/d.numLatticeSites)
      blockcoord[2] = math.floor((i.coordinates[2])/d.numLatticeSites)
      blockID=int((blockcoord[0]*d.blockSize[1]+blockcoord[1])*d.blockSize[2]+blockcoord[2])

      siteID=int((i.coordinates[0]%d.numLatticeSites)*d.numLatticeSites+(i.coordinates[1]%d.numLatticeSites))*d.numLatticeSites+(i.coordinates[2]%d.numLatticeSites)
      cpus = []

      for key in entries:
        core = key[0]
        block = key[1]
        if (block==blockID):
          cpus.append(core)

      p.pack_uint(d.dictIdx[i.vertexID])
      p.pack_uint(i.coreNum)
      p.pack_uint(i.coordinates[0])
      p.pack_uint(i.coordinates[1])
      p.pack_uint(i.coordinates[2])
      p.pack_uint(i.siteType)
      p.pack_uint(len(i.biEdges))
 
      for j in i.biEdges:
          x = d.dictIdx.get(j,-1)
          if (x == -1): continue
          p.pack_uint(d.dictIdx[j])
  
      p.pack_uint(blockID) 
      p.pack_uint(d.blockFluidSiteCounts[blockID])
      p.pack_uint(siteID)
   
      p.pack_uint(len(cpus))

      for j in cpus:
        p.pack_uint(j)

    f.write(p.get_buffer())       

def writeForCBinPlus4000(d, path):
    entries = []
    f = open(path,'wb')
    blockcoord = [0.0,0.0,0.0]

    p = xdrlib.Packer()

    for i in d.vertices:
      # Read halo- blocks on domain edges
      if (i.coordinates[0]%d.numLatticeSites == 0): #Lower x-edge
        xRange = [-1,0]
      elif (i.coordinates[0]%d.numLatticeSites == 7): #Upper x-edge
        xRange = [0,1]
      else:
        xRange = [0]

      if (i.coordinates[1]%d.numLatticeSites == 0): #Lower x-edge
        yRange = [-1,0]
      elif (i.coordinates[1]%d.numLatticeSites == 7): #Upper x-edge
        yRange = [0,1]
      else:
        yRange = [0]

      if (i.coordinates[2]%d.numLatticeSites == 0): #Lower x-edge
        zRange = [-1,0]
      elif (i.coordinates[2]%d.numLatticeSites == 7): #Upper x-edge
        zRange = [0,1]
      else:
        zRange = [0]


      for dx in xRange:
        for dy in yRange:
          for dz in zRange:
            blockcoord[0] = math.floor((i.coordinates[0]+dx)/d.numLatticeSites)
            blockcoord[1] = math.floor((i.coordinates[1]+dy)/d.numLatticeSites)
            blockcoord[2] = math.floor((i.coordinates[2]+dz)/d.numLatticeSites)
            blockid=int((blockcoord[0]*d.blockSize[1]+blockcoord[1])*d.blockSize[2]+blockcoord[2])

            key = (i.coreNum,blockid)
            if key not in entries:
                entries.append(key)
	        print [i.coreNum,blockid,d.blockOffset[blockid],d.blockDataLength[blockid],d.blockUncompressedDataLength[blockid]]
                p.pack_uint(i.coreNum)
	        p.pack_uint(blockid)
		p.pack_uint(d.blockOffset[blockid])
		p.pack_uint(d.blockDataLength[blockid])
		p.pack_uint(d.blockUncompressedDataLength[blockid])
	  	p.pack_uint(d.blockFluidSiteCounts[blockid])
     
    f.write(p.get_buffer())


def readForCBin(path):
    """
    Read a file from the binary format, provided a path
    Core ID - Coordinates - SiteTypeNum - Edges - -1
    """
    f = open(path,'r')

    d = domain.Domain() 

    # 1. Total site aocunt
    # x = struct.pack('i',totalSites)

    size_int = 4
    size_double = 8

    x = f.read(size_int)
    len_vertices = struct.unpack('i',x)[0]
    #print "Site count = ", len_vertices
    d.numLatticeSites = 8 

    len_biEdges = []
    # 2. neighbour counts for each Vertex
    for i in xrange(0, len_vertices):
        v = domain.VertexGmy()
        # x = struct.pack('i',i.vertexID)
        # x = struct.pack('i',v.dictIdx[i.vertexID])
        v.vertexID = struct.unpack('i',f.read(size_int))[0]
        # x = struct.pack('i',len(i.biEdges))
        len_biEdges.append(struct.unpack('i',f.read(size_int))[0])
        d.vertices.append(v)
        d.dictIdx[i] = v.vertexID
        # print "read: ", v.vertexID, len_biEdges

 
    vertex_iter = 0

    # 3. Vertex data 
    for i in d.vertices:
        # x = struct.pack('i',i.vertexID)

        #TODO: fix!
        i.vertexID = struct.unpack('i',f.read(size_int))[0]
        # x = struct.pack('i',v.dictIdx[i.vertexID])

        i.coreNum = struct.unpack('i',f.read(size_int))[0]
        # x = struct.pack('i',i.coreNum)

        coord0 = struct.unpack('i',f.read(size_int))[0]
        # x = struct.pack('i',i.coordinates[0])
        coord1 = struct.unpack('i',f.read(size_int))[0]
        # x = struct.pack('i',i.coordinates[1])
        coord2 = struct.unpack('i',f.read(size_int))[0]
        # x = struct.pack('i',i.coordinates[2])
        i.coordinates = [coord0, coord1, coord2]

        i.siteType = struct.unpack('i',f.read(size_int))[0]
        # x = struct.pack('i',i.siteType)
        
        i.biEdges = []
        for j in xrange(0, len_biEdges[vertex_iter]):
            #TODO: check!
            i.biEdges.append(struct.unpack('i',f.read(size_int))[0])
            # x = struct.pack('i',v.dictIdx[j])
    
        print "read: ", i.vertexID, i.coreNum, i.coordinates, i.siteType
        vertex_iter += 1

    d.updateProperties()

    return d

def fix_xadj(v):
    """
    Converts the xadj list to a format that is in line with the CSR standard.
    """
    for c in xrange(len(v.xadj)):
        # Prepend a 0 to each array (each core).
        v.xadj[c] = [0] + v.xadj[c]

        # Make the xadj counting cumulative.
        count = 0
        for i in xrange(len(v.xadj[c])):
            v.xadj[c][i] += count             
            count = v.xadj[c][i]        
            
    return v

def AddVertex(c, vertex):
    """
    Adds a vertex to the xadj list (part of the CSR format).
    """
    core = vertex.coreNum
    
    #if (isinstance(i,domain.VertexGmy)):
    #    v.xadj[c]     += [len(i.biEdges)] 
    #else:
    c.xadj[core]     += [len(vertex.edges)]
    
    # TODO: Ignoring the one-directional edges for now, but this could become an issue 
    # when we write GMYs!
    
    # we only add neighbour counters to xadj.
    # xadj will need to be fixed to make it ascending and to prepend a 0.
    
    c.vids[core] += [vertex.vertexID]
    c.vcoords[core] += vertex.coordinates

    local_adjncy = []    
    # GMY edges are tuples of (destId, edgetype). We are interested in destId
    # here
    # edgetype could matter in the future, if we wish to incorporate edge
    # weighting.
    for j in vertex.edges: 
        local_adjncy.append(j)
        
    c.adjncy[core] += local_adjncy
    c.wgts[core]     += [vertex.siteType]
        
    return c

def ConvertToCSR(d):
    """ Converts the Domain Object to a CSR structure
    """
    
    #self.vertices = []             #List of all vertices
    #self.numGraphs = 0             #Number of subgraphs, i.e. number of processors
    #self.numSites = []             #Number of sites (Number of vertices per subgraph)
    #self.neighCount = []         #Neighbor Count (Number of subgraphs linked to this subgraph)
    #self.commSurf = []             #Communication Surface (Number of edges crossing subgraph boundary) 
    #self.numEdges = []             #Number Of Edges (Number of edges within a subgraph)
    #self.totalWeight = []        #Total Weight according to site Type (Gmy)
    #self.dictIdx = {}                #Map vertexID to its index in the vertices list
    
    c = csr.CSR()    
    
    c.vtxdist = [0]
    
    c.xadj = [[] for x in xrange(d.numGraphs)]    
    c.adjncy = [[] for x in xrange(d.numGraphs)]
    c.wgts = [[] for x in xrange(d.numGraphs)] 
    c.vids = [[] for x in xrange(d.numGraphs)]
    c.vcoords = [[] for x in xrange(d.numGraphs)]
    
    for core in xrange(0, d.numGraphs):
        c.vtxdist += d.numSites[core]
        
    #local2global_site_id_mapping = [site_id on c]
    for i in d.vertices:
        c = AddVertex(c, i)

    c = fix_xadj(c)
                     
    return c

def ConvertFromCSR(c):
    part_loc = c
    #PPStee output arrays (each core has a simple array of length <num sites> and each entry assigns a local site to a target core, 
    # i.e. part_loc[i] = core id where the local site i was assigned to)
#    for c in xrange(0,core_number):
#    for i, target in part_loc
#        global_site_id = local2global_site_id_mapping
#        target_core = target
        #write to binary file: 
#        x = struct.pack('ii', global_site_id, target_core, other_site_data(current_global_site_id)
#        f.write(x)

    
def WriteCSRFile(CSRclass, fname, XDR_format="no"):
    """ Writes a CSR class to file
    """
    c = CSRclass

    num_cores = len(c.vtxdist) - 1
    
    f = open(fname, "w")    
    
    if XDR_format == "no":
             
        x = struct.pack('i',num_cores)
        f.write(x)
    
        # Pack vtxdist
        x = struct.pack('i',0)
        f.write(x)
        for i in xrange(0, num_cores):
            x = struct.pack('i', c.vtxdist[i])
            f.write(x)
 
            # For each core, pack vids, vcoords, vwgts, xadj and adjncy.
        for i in xrange(0, num_cores):
            length_of_vertexdata = c.vtxdist[i+1] -    c.vtxdist[i]
            #print length_of_vertexdata, len(c.vids[i])
            
            for j in xrange(0, length_of_vertexdata):
                x = struct.pack('i', c.vids[i][j])
                f.write(x)
    
            for j in xrange(0, length_of_vertexdata*3):
                x = struct.pack('i', c.vcoords[i][j])
                f.write(x)
 
            for j in xrange(0, length_of_vertexdata):
                x = struct.pack('i', c.wgts[i][j])
                f.write(x)
                         
            x = struct.pack('i', 0)
            f.write(x)
            for j in xrange(0, length_of_vertexdata):
                x = struct.pack('i', c.xadj[i][j])
                f.write(x)

            for j in xrange(0, len(c.adjncy[i])):
                #print i, len(c.adjncy), j, len(c.adjncy[i])
                
                #print c.adjncy[i]                                
                
                x = struct.pack('i', c.adjncy[i][j])
                f.write(x)                
        
    else:
        p = xdrlib.Packer()    
    
        p.pack_uint(num_cores)
    
        # Pack vtxdist
        p.pack_uint(0)
        for i in xrange(0, num_cores):
            p.pack_uint(c.vtxdist[i])
 
        # For each core, pack vids, vcoords, vwgts, xadj and adjncy.
        for i in xrange(0, num_cores):
            length_of_vertexdata = c.vtxdist[i+1] -    c.vtxdist[i]
            #print length_of_vertexdata, len(c.vids[i])
            for j in xrange(0, length_of_vertexdata):
                p.pack_uint(c.vids[i][j]) 
                    
            for j in xrange(0, length_of_vertexdata*3):
                p.pack_uint(c.vcoords[i][j]) 
        
            for j in xrange(0, length_of_vertexdata):
                p.pack_uint(c.vwgts[i][j])             
            
            p.pack_uint(0)
            for j in xrange(0, length_of_vertexdata):
                p.pack_uint(c.xadj[i][j])            

            for j in xrange(0, length_of_vertexdata):
                p.pack_uint(c.adjncy[i][j])            
    
        #Write the XDR buffer
        f.write(p.get_buffer())
        
    return c

def ReadCSRFile(fname, XDR_format="no"):
    f = open(fname,'r')
    c = csr.CSR()    
    
    
    if XDR_format == "yes":
        compressed = f.read()
        #uncompressed = zlib.decompress(compressed) #TODO: Add compression if we want that...
        u = xdrlib.Unpacker(compressed)
        
        """
        CSR format:
        All values are uints.
        
        num_cores 
        vtxdist: distribution of vertices over processes (has length num_cores+1)
     
        then, for each core:
        length_of_vertexdata (can be obtained from vtxdist)
        vids: vertexID list for each proc                             LEN:length_of_vertexdata
        vcoords: vertex coordinate list for each proc    LEN:length_of_vertexdata*3
        vwgts: vertex weights                                                    LEN:length_of_vertexdata
        xadj: neighbour list offsets                                     LEN:length_of_vertexdata+1
        adjncy: adjacency list (uses global vars)            LEN:xadj[-1]
     
        """
        
        num_cores = u.unpack_uint()
        vtxdist = np.array(([0]))
        for i in xrange(0, num_cores):
            c.vtxdist = np.append(vtxdist,u.unpack_uint())
     
        for i in xrange(0, num_cores):
            length_of_vertexdata = vtxdist[i+1] -    vtxdist[i]
            vids = np.array(([]))
            for j in xrange(0, length_of_vertexdata):
                vids = np.append(vids, u.unpack_uint())
            c.vids = c.vids + [vids]
                
            vcoords = np.array(([]))
            for j in xrange(0, length_of_vertexdata*3):
                vcoords = np.append(vcoords, u.unpack_uint())
            c.vcoords = c.vcoords + [vcoords]
            
            vwgts = np.array(([]))
            for j in xrange(0, length_of_vertexdata):
                vwgts = np.append(vwgts, u.unpack_uint())
            c.vwgts = c.vwgts + [vwgts]
            
            xadj = np.array(([0]))
            for j in xrange(0, length_of_vertexdata):
                xadj = np.append(xadj, u.unpack_uint())
            c.xadj = c.xadj + [xadj]
            
            adjncy = np.array(([]))
            for j in xrange(0, length_of_vertexdata):
                adjncy = np.append(adjncy, u.unpack_uint())
            c.adjncy = c.adjncy + [adjncy]

    return c
