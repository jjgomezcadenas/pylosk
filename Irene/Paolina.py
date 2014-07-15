"""
Calls Paolina
"""
from ROOT import *
import numpy as np

# Load paolina library
gSystem.Load("libpaolina")

from KLog import *
#create logger
lgx =logging.getLogger("Paolina")
lgx.setLevel(logging.DEBUG)
lgx.addHandler(ch)
debug = Debug.verbose.value

class PVoxel(object):
    """
    A wrapper of Paolina Voxels in python 
    JJ, Summer, 2014
    """

    def __init__(self,pvoxel):
        """
        Constructor takes a pvoxel : The Paolina C++ voxel object
        """

        pos = pvoxel.Pos()
        apos = pvoxel.AbsPos()
        size = pvoxel.Size()

        xpos =[]
        xapos =[]
        xsize=[]

        for i in range(0,3):
            xpos[i]=pos[i]
            xapos[i]=apos[i]
            xsize[i]=size[i]

        self.ipos = np.array(xpos)
        self.pos = np.array(xapos)
        self.size = np.array(size)
        self.edep = pvoxel.Edep()
    def Position(self):
        """
        Returns the position in mm
        """
        return self.pos

    def Edep(self):
        """
        Returns the energy in MeV
        """
        return self.edep

    def Size(self):
        """
        Returns the size in mm
        """
        return self.size
 
    


class Paolina(object):
    """
    Acess to Paolina 
    JJ, Summer, 2014
    """

    def __init__(self,voxelXYZ,blobRadius=20):
    	"""
    	Constructor takes:
        ihits: a vector of hits (eventually smeared)
        voxelXYZ = (x,y,z) size of voxel (in mm) 
        blobRadius: radius of blob  (in mm)

    	"""

        self.Rblob = blobRadius
        self.vXYZ = voxelXYZ
        pvoxel_size = std.vector("double")()
    

        pvoxel_size.push_back(voxelXYZ[0]) #voxel size
        pvoxel_size.push_back(voxelXYZ[1])
        pvoxel_size.push_back(voxelXYZ[2])

        left_range = std.vector("double")()
        right_range = std.vector("double")()
        left_range.push_back(-500.)
        left_range.push_back(-500.)
        left_range.push_back(-500.)
        right_range.push_back(500.)
        right_range.push_back(500.)
        right_range.push_back(500.)

        self.pTB = paolina.TrackBuilder() # Creates a Paolina track builder
        self.pVB =paolina.VoxelBuilder(pvoxel_size, left_range,right_range) #voxel builder
        self.pBB = paolina.BlobBuilder(blobRadius) #Blob Builder


    def Voxels(self,ievt):
        """
        Returns voxels: takes a handel to irene event
        Notice that ihit is a C++ object of type
        std::vector<std::vector<double> >
        Paolina returns objects:
        std::vector<paolina::Voxel*> voxels
        """
        ihits = std.vector("std::vector<double>")()
        ihits = ievt.FillHitVector("ACTIVE") #get all hits in event

        for i in range (0,ihits.size()):
            hit = ihits[i]
            print "hit x=%7.1f y=%7.1f z=%7.1f e=%7.1e "%(hit[0],hit[1],hit[2],hit[3])

        cond_pause(Debug.verbose.value)

        voxels = self.pVB.GetVoxels(ihits)
        ListVox =[]

        for i in range(0,voxels.size()):
            ListVox.append(PVoxel(voxels.at(i)))

        return ListVox

  
  