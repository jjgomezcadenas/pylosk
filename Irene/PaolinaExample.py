"""
An example of how to call Paolina: 
"""


from IEventReader import IEventReader
from IParticle import IParticle
from Paolina import Paolina, PVoxel
from TDrawHits import TDrawHits
from MPLDrawHits import MPLDrawHits
import IParam as ip
import sys,getopt
import numpy as np

from KLog import *
#create logger
lgx =logging.getLogger("PaolinaExample")
lgx.setLevel(logging.DEBUG)
lgx.addHandler(ch)
debug = Debug.verbose.value
DrawRoot = False
DrawMPL = False

def main(argv):
    
    pathToFile,nEvents = GetArguments(argv)

    #reader---
    eventReader = IEventReader(pathToFile)
    nRun = min(eventReader.NumberOfEvents(),nEvents)
    
    #--Logging
    s = "Starting Event Reader: number of events requested ={0}"
    s+=" number of events in file ={1}:"
    s+=" number of events to run ={2}"
    lgx.info(s.format(nEvents,eventReader.NumberOfEvents(),nRun))
    cond_pause(Debug.quiet.value)
    #--

    #Init Paolina
    voxelXYZ=[1.,1.,1.] # in mm
    blobRadius=20 # in mm
    plna = Paolina(voxelXYZ,blobRadius)

  #------Loop
    for event in range(0,nRun):
    #---------

        #--Logging
        s="Reading event number {0}"
        lgx.debug(s.format(event))
        cond_pause(Debug.verbose.value)
        #--

        #read event
        ievt = eventReader.ReadEvent(event)

        #--Logging
        s="Number of Bytes read = {0} "
        s+="Total number of Bytes read = {1}"
        lgx.debug(s.format(eventReader.NumberOfBytesRead(),
                        eventReader.TotalNumberOfBytesRead()))
        cond_pause(Debug.verbose.value)
        #--

        #Get MC particle
        smearVector =np.array([ip.sigma_x, ip.sigma_y, 0.]) 
        betaMax = IParticle(ievt,smearVector,ip.sample)
        trueHits = betaMax.TrueHits()
        smearHits = betaMax.SmearedHits()

        #--Logging
        s ="betaMax ={0}"
        lgx.debug(s.format(betaMax))
        cond_pause(Debug.verbose.value)
        if debug >= Debug.verbose.value:
            DrawHits(trueHits,smearHits)
        #----

        #Call Paolina
        voxels = plna.Voxels(ievt)
        

       
    

        

def GetArguments(argv):
    inputFile = ''
    inputDir = ''
    outputDir = ''
    n_events=''
    try:      
        opts, args = getopt.getopt(argv,"hi:d:o:n:",["ifile=","idir=","odir=","events="])
    except getopt.GetoptError:
        print 'IMain -i <inputFile> -d <inputDir> -n <max number of events> -g <generator>'
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print 'pylosk -i <inputFile> -d <inputDir> -o <outputDir> -n <max number of events>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputFile = arg
        elif opt in ("-d", "--idir"):
            inputDir = arg
        elif opt in ("-n", "--events"):
            n_events = arg
        elif opt in ("-o", "--odir"):
            outputDir = arg
            
    if inputFile=="":
        inputFile=ip.inputFile
    if inputDir=="":
        inputDir=ip.inputDir
    if outputDir=="":
        outputDir=ip.inputDir

    if n_events=="":
        n_events=ip.n_events
    else:
        n_events=int(n_events)

    pathToFile = inputDir+'/'+inputFile

    s = "Reading file ={0}"
    lgx.info(s.format(pathToFile))
    
    return (pathToFile,n_events)

def DrawHits(trueHits,smearHits):
    """
    Given a KMCParticle, draw hits, using root or Matplotlib
    """

    print "drawing"
    
    scale =.6*ip.Pr
    lowerBin =np.array([-700.,-700.,-700.])/scale
    upperBin =np.array([700.,700.,700.])/scale

    if DrawRoot == True:
        drawHits = TDrawHits((lowerBin,upperBin),trueHits,smearHits)
        drawHits.DrawMeasurements(draw='2D')
        drawHits.DrawMeasurements(draw='3D')

        drawHits.DrawHits(draw='2D')
        drawHits.DrawHits(draw='3D')

        drawHits.DrawAll(draw='2D')
        
    elif DrawMPL == True:
        drawHits = MPLDrawHits((lowerBin,upperBin),trueHits,smearHits)
        # drawHits.DrawMeasurements(draw='2D')
        # drawHits.DrawMeasurements(draw='3D')

        drawHits.DrawHits(draw='2D')
        drawHits.DrawHits(draw='3D')

        #drawHits.DrawAll(draw='2D')
        

 

if __name__ == '__main__':
    main(sys.argv[1:])


