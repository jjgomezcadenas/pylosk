"""
Main driver: a file with event and calls the user-supplied method: 
"""


from IEventReader import IEventReader
from IParticle import IParticle
from KTrackFitter import KTrackFitter
from KFWolinFilter import KFWolinFilter
from TDrawHits import TDrawHits
from MPLDrawHits import MPLDrawHits
import IParam as ip
import sys,getopt
import numpy as np

from KLog import *
#create logger
lgx =logging.getLogger("IMain")
lgx.setLevel(logging.DEBUG)
lgx.addHandler(ch)
debug = Debug.verbose.value
DrawRoot = False
DrawMPL = True

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

    #Kalman Filter Fitter
    slkf = KFWolinFilter("KFWolin",P0=2.9) #KF Fitter 
    
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
        

        # Create a KalmanFilter to be used to do the fitting
        lgx.debug("-- Creating KFWolinFilter...")
        kfilter = KFWolinFilter("KFWolinFilter")

        #Define a track fitter
        lgx.debug("-- Creating KFTrackFitter...")
        tfitter  = KTrackFitter(kfilter,smearHits,smearVector,
                                 betaMax,chi2Limit=ip.Chi2Lim,Pressure=ip.Pr)
    
        
        # Perform the fit.
        lgx.debug("-- Performing fit...")
        tfitter.Fit()
        cond_pause(Debug.verbose.value)

        segments = tfitter.Segments

        lgx.debug("-- Fitter found {0} segments".format(segments))
    
        for seg in segments:
        
            # Get the mean, min, and max chi2 values for each segment.
            #  Set to -1 if the segment is not >= 3 points.
            
            mean_chi2 = -1.; min_chi2 = -1.; max_chi2 = -1.;
            
            if(len(seg.seg_k) >= 2):
                mean_chi2 = np.mean(seg.seg_fchisq[1:])
                min_chi2 = min(seg.seg_fchisq[1:])
                max_chi2 = max(seg.seg_fchisq[1:])
        
            # Print the segment file line.
            s="segment id = {0} length ={1}"
            s+=" mean chi2 ={2} min chi2 = {3} max chi2 ={4}\n"
            lgx.info(s.format(seg.seg_id,
            len(seg.seg_k),mean_chi2,min_chi2,max_chi2));
        
       
    

        

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


