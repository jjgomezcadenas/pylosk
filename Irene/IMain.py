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

from KFBase import KFVector
from KLog import *
#create logger
lgx =logging.getLogger("IMain")
lgx.setLevel(logging.WARN)
lgx.addHandler(ch)
#debug = Debug.verbose.value
debug = Debug.mute.value
DrawRoot = False
DrawMPL = True

# Temporary path variables for fit writing
import os
fit_outdir="/data4/NEXT/users/jrenner/kalmanfilter/out/ifit"

def main(argv):
    
    pathToFile,sevt,eevt,itrk_name,bbevt,rev_trk = GetArguments(argv)

    print "-- Got args.  pathToFile = {0}, sevt = {1}, eevt = {2}, itrk_name = {3}, bbevt = {4}, rev_trk = {5}".format(pathToFile,sevt,eevt,itrk_name,bbevt,rev_trk);

    if(rev_trk):
        print "\n\n-- WORKING ON REVERSED TRACKS --\n\n"
        #fnb_trk = "{0}/rev/".format(trk_outdir)
        fnb_fit = "{0}/{1}/rev/".format(fit_outdir,itrk_name)
    else:
        fnb_fit = "{0}/{1}".format(fit_outdir,itrk_name)

    if(not os.path.isdir(fit_outdir)):
        try:
            os.mkdir(fit_outdir);
        except OSError:
            print "Error creating {0}".format(fit_outdir);
    if(not os.path.isdir("{0}/{1}".format(fit_outdir,itrk_name))):
        print "Directory {0}/{1} does not exist, creating...".format(fit_outdir,itrk_name);
        try:
            os.mkdir("{0}/{1}".format(fit_outdir,itrk_name));
        except OSError:
            print "Error creating {0}/{1}".format(fit_outdir,itrk_name);
    if(not os.path.isdir("{0}/{1}/rev".format(fit_outdir,itrk_name))):
        print "Directory {0}/{1}/rev does not exist, creating...".format(fit_outdir,itrk_name);
        try:
            os.mkdir("{0}/{1}/rev".format(fit_outdir,itrk_name));
        except:
            print "Error creating {0}/{1}/rev".format(fit_outdir,itrk_name);

    #reader---
    eventReader = IEventReader(pathToFile)
    #nRun = min(eventReader.NumberOfEvents(),nEvents)
    
    #--Logging
    #s = "Starting Event Reader: number of events requested ={0}"
    #s+=" number of events in file ={1}:"
    #s+=" number of events to run ={2}"
    #lgx.info(s.format(nEvents,eventReader.NumberOfEvents(),nRun))
    #cond_pause(Debug.quiet.value)
    #--

    #Kalman Filter Fitter
    #slkf = KFWolinFilter("KFWolin",P0=2.9) #KF Fitter 
    
    #------Loop: cover events sevt to (eevt-1)
    for event in range(sevt,eevt):
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
        smearVector =np.array([ip.sigma_x, ip.sigma_y, 0]) 
        betaMax = IParticle(ievt,smearVector,ip.sample)
        
        # Get the other electron if this is a double-beta event.
        if(bbevt):
            betaSecond = IParticle(ievt,smearVector,ip.sample,2447.,2)

            # Assign the list of smeared hits appropriately.            
            if(not rev_trk):
                trueHits = betaMax.TrueHits(); trueHits.reverse()
                secondHits = betaSecond.TrueHits()
                smearHits = betaMax.SmearedHits(); smearHits.reverse()
                sSmearHits = betaSecond.SmearedHits()
            else:
                trueHits = betaSecond.TrueHits(); trueHits.reverse()
                secondHits = betaMax.TrueHits()
                smearHits = betaSecond.SmearedHits(); smearHits.reverse()
                sSmearHits = betaMax.SmearedHits()
            
            # Construct the list of true hits from both tracks.
            for th in secondHits:
                trueHits.append(th)
                
            # Construct the list of smeared hits from both tracks.
            for sh in sSmearHits:
                smearHits.append(sh)
        else:
            trueHits = betaMax.TrueHits()
            smearHits = betaMax.SmearedHits(rev_trk)                        

#        for th in trueHits:
#            print "{0} {1} {2} 0".format(th[0],th[1],th[2]);
#        for sh in secondHits:
#            print "{0} {1} {2} 1".format(sh[0],sh[1],sh[2]);
        
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
        
        # Temporary fit result writing
        f_ftrk = open("{0}/fit_{1}_{2}.dat".format(fnb_fit,itrk_name,event),"w")
        f_fseg = open("{0}/seg_{1}_{2}.dat".format(fnb_fit,itrk_name,event),"w")
        f_ftrk.write("# segID k x0 y0 z0 p1p p2p p3p p4p chi2p p1f p2f p3f p4f chi2f cfxy cftxy\n")
        f_fseg.write("# segID nPts chi2avg chi2min chi2max\n")

        for seg in segments:
        
            # Get the mean, min, and max chi2 values for each segment.
            #  Set to -1 if the segment is not >= 3 points.
            mean_chi2 = -1.; min_chi2 = -1.; max_chi2 = -1.;
            if(len(seg.seg_k) >= 2):
                mean_chi2 = np.mean(seg.seg_fchisq[1:]);
                min_chi2 = min(seg.seg_fchisq[1:]);
                max_chi2 = max(seg.seg_fchisq[1:]);
            
            # Print the segment file line.
            f_fseg.write("{0} {1} {2} {3} {4}\n".format(seg.seg_id,len(seg.seg_k),mean_chi2,min_chi2,max_chi2));
            
            # Print the track file lines for each segment point: include the smeared hits.
            for k,x0,y0,z0,p1p,p2p,p3p,p4p,chi2p,p1f,p2f,p3f,p4f,chi2f,cfxy,cftxy in zip(seg.seg_k,seg.seg_x0,seg.seg_y0,seg.seg_z0,seg.seg_p1p,seg.seg_p2p,seg.seg_p3p,seg.seg_p4p,seg.seg_pchisq,seg.seg_p1f,seg.seg_p2f,seg.seg_p3f,seg.seg_p4f,seg.seg_fchisq,seg.seg_cfxy,seg.seg_cftxy):
                if(isinstance(chi2f,KFVector)):
                    chi2 = chi2f[2];
                else:
                    chi2 = chi2f;
                f_ftrk.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12} {13} {14} {15} {16}\n".format(seg.seg_id,k,x0,y0,z0,p1p,p2p,p3p,p4p,chi2p,p1f,p2f,p3f,p4f,chi2f,cfxy,cftxy));
    

        

def GetArguments(argv):
    inputFile = ''
    inputDir = ''
    outputDir = ''
    sevt=''
    eevt=''
    bbevt=False
    rev_trk=False
    itrk_name='' 
    try:      
        opts, args = getopt.getopt(argv,"hi:d:s:e:u:g:r:",["ifile=","idir=","sevt=","eevt","rname=","gen=","rev="])
    except getopt.GetoptError:
        print 'IMain -i <inputFile> -d <inputDir> -s <event start> -e <event end> -u <run_name> -g <0=sel,1=bb> -r <0=forward,1=reverse>'
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            print 'IMain -i <inputFile> -d <inputDir> -s <event start> -e <event end> -u <run_name> -g <0=sel,1=bb> -r <0=forward,1=reverse>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputFile = arg
        elif opt in ("-d", "--idir"):
            inputDir = arg
        elif opt in ("-s", "--sevt"):
            sevt = arg
        elif opt in ("-e", "--eevt"):
            eevt = arg
        elif opt in ("-o", "--odir"):
            outputDir = arg
        elif opt in ("-u", "--rname"):
            itrk_name = arg
        elif opt in ("-g", "--gen"):
            if(int(arg) > 0): bbevt = True
            else: bbevt = False
        elif opt in ("-r", "--rev"):
            if(int(arg) > 0): rev_trk = True
            else: rev_trk = False 
            
    if inputFile=="":
        inputFile=ip.inputFile
    if inputDir=="":
        inputDir=ip.inputDir
    if outputDir=="":
        outputDir=ip.inputDir
    if itrk_name=="":
        itrk_name = "sel"

    if (sevt=="" or eevt==""):
        sevt=0
        eevt=ip.n_events
    else:
        sevt=int(sevt)
        eevt=int(eevt)

    pathToFile = inputDir+'/'+inputFile

    return (pathToFile,sevt,eevt,itrk_name,bbevt,rev_trk)

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


