"""
kftrackfit.py

Attempts a Kalman filter fit on a track using KFTrackFitter

Fit output format:
    
    fit file: k p1p p2p p3p p4p chi2p p1f p2f p3f p4f chi2f
    segment file:
        
    All distances are in cm

"""
import sys
import numpy as np
import scipy.integrate as integrate
import random as rd
import os
from math import *
from trackdefs import *
from scipy.interpolate import interp1d

from ToyParticle import ToyParticle
from KTrackFitter import KTrackFitter
from KFBase import KFVector
from KFWolinFilter import KFWolinFilter

from KLog import *
#create logger
lgx =logging.getLogger("ktrackfit")
lgx.setLevel(logging.INFO)
lgx.addHandler(ch)
debug = Debug.info.value

# File names.
if(rev_trk):
    print "\n\n-- WORKING ON REVERSED TRACKS --\n\n"
    #fnb_trk = "{0}/rev/".format(trk_outdir)
    fnb_trk = "{0}/{1}/".format(trk_outdir,trk_name)
    fnb_fit = "{0}/{1}/rev/".format(fit_outdir,trk_name)
else: 
    fnb_trk = "{0}/{1}".format(trk_outdir,trk_name)
    fnb_fit = "{0}/{1}".format(fit_outdir,trk_name)

if(not os.path.isdir("{0}/{1}".format(trk_outdir,trk_name))):
    print "ERROR: no tracks available in {0}/{1}".format(trk_outdir,trk_name);
    sys.exit();
if(not os.path.isdir(fit_outdir)): os.mkdir(fit_outdir);
if(not os.path.isdir("{0}/{1}".format(fit_outdir,trk_name))): os.mkdir("{0}/{1}".format(fit_outdir,trk_name));
if(not os.path.isdir("{0}/{1}/rev".format(fit_outdir,trk_name))): os.mkdir("{0}/{1}/rev".format(fit_outdir,trk_name));

# Create num_tracks tracks.
for ntrk in range(num_tracks):

    logging.info("\n\n-- Track {0} --\n\n".format(ntrk))

    # Create a ToyParticle for this track.
    logging.debug("-- Creating ToyParticle...")
    tfile = "{0}/{1}_{2}.dat".format(fnb_trk,trk_name,ntrk)
    tpart = ToyParticle(tfile,rev_trk,np.array([0.5,0.5,0.5]),0)
    
    # Create a KalmanFilter to be used to do the fitting
    logging.debug("-- Creating KFWolinFilter...")
    kfilter = KFWolinFilter("KFWolinFilter");
    
    # Set up a KFTrackFitter to fit this track.
    logging.debug("-- Creating KTrackFitter...")
    tfitter = KTrackFitter(kfilter,tpart.SmearedHits(rev_trk),np.array([0.5,0.5]),tpart,chi2_lim,10.)
    
    # Perform the fit.
    logging.debug("-- Performing fit...")
    tfitter.Fit()
        
    # Write the fit files.
    f_ftrk = open("{0}/fit_{1}_{2}.dat".format(fnb_fit,trk_name,ntrk),"w")
    f_fseg = open("{0}/seg_{1}_{2}.dat".format(fnb_fit,trk_name,ntrk),"w")
    f_ftrk.write("# segID k x0 y0 z0 p1p p2p p3p p4p chi2p p1f p2f p3f p4f chi2f cfxy cftxy\n")
    f_fseg.write("# segID nPts chi2avg chi2min chi2max\n")
    
    segments = tfitter.Segments
    
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

    # Close the files.
    f_ftrk.close();
    f_fseg.close();
