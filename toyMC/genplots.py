"""
genplots.py

Plots the relevant statistical quantities for the forward and reverse tracks.

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from math import *
from trackdefs import *

from abc import ABCMeta, abstractmethod
import logging 

# Cut parameters
mcut = 0.05;
bcut = -0.5;

smcut = -1.4;
sbcut = 0.9;

# Ensure the correct directories exist.
if(not os.path.isdir(fit_outdir)): 
    print "ERROR: fit directory {0} not available".format(fit_outdir);
    sys.exit();
if(not os.path.isdir("{0}/{1}/rev".format(fit_outdir,trk_name))): 
    print "ERROR: fit directory {0}/rev is not available".format(fit_outdir);
    sys.exit();
if(not os.path.isdir("{0}/{1}".format(fit_outdir,trk_name))):
    print "ERROR: fit directory {0}/{1} is not available".format(fit_outdir,trk_name);
    sys.exit();
if(not os.path.isdir(plt_outdir)): 
    print "ERROR: plot directory {0} not available".format(plt_outdir);
    sys.exit();
if(not os.path.isdir("{0}/{1}".format(plt_outdir,trk_name))):
    os.mkdir("{0}/{1}".format(plt_outdir,trk_name));
    print "Created plot directory {0}/{1}".format(plt_outdir,trk_name);

plt_base = "{0}/{1}".format(plt_outdir,trk_name);

# Plot options.
plt_show = False;
plt_print = True;

# Extract data from the segments files for num_tracks tracks.
splot_ktot = []; splot_rktot = [];
splot_chi2avg = []; splot_rchi2avg = [];
splot_sfac = []; splot_rsfac = [];
splot_schi2_ff = []; splot_schi2_fl = [];
splot_schi2_rf = []; splot_schi2_rl = [];
splot_ptps_ff = []; splot_ptps_fl = [];
splot_ptps_rf = []; splot_ptps_rl = [];
splot_scfxy_ff = []; splot_scfxy_fl = [];
splot_scfxy_rf = []; splot_scfxy_rl = [];
splot_cfxy_min = []; splot_cfxy_konmin = [];
splot_rcfxy_min = []; splot_rcfxy_konmin = [];

prof_fkon = []; prof_fchi2 = [];
prof_rkon = []; prof_rchi2 = [];
prof_fckon = []; prof_fcfxy = [];
prof_rckon = []; prof_rcfxy = [];
for ntrk in range(num_tracks):
    
    logging.debug("-- Segment analysis for track {0}\n".format(ntrk));
        
    # Read in the forward fit.
    # segID k p1p p2p p3p p4p chi2p p1f p2f p3f p4f chi2f
    fittbl = np.loadtxt("{0}/{1}/fit_{2}_{3}.dat".format(fit_outdir,trk_name,trk_name,ntrk));
    fit_seg = fittbl[:,0];
    fit_k = fittbl[:,1];
    fit_x0 = fittbl[:,2];
    fit_y0 = fittbl[:,3];
    fit_z0 = fittbl[:,4];
    fit_p1p = fittbl[:,5];
    fit_p2p = fittbl[:,6];
    fit_p3p = fittbl[:,7];
    fit_p4p = fittbl[:,8];
    fit_chi2p = fittbl[:,9];
    fit_p1f = fittbl[:,10];
    fit_p2f = fittbl[:,11];
    fit_p3f = fittbl[:,12];
    fit_p4f = fittbl[:,13];
    fit_chi2f = fittbl[:,14];
    fit_cfxy = fittbl[:,15];
    fit_cftxy = fittbl[:,16];   
 
    # Read in the reverse fit.
    # segID k p1p p2p p3p p4p chi2p p1f p2f p3f p4f chi2f
    rfittbl = np.loadtxt("{0}/{1}/rev/fit_{2}_{3}.dat".format(fit_outdir,trk_name,trk_name,ntrk));
    rfit_seg = rfittbl[:,0];
    rfit_k = rfittbl[:,1];
    rfit_x0 = rfittbl[:,2];
    rfit_y0 = rfittbl[:,3];
    rfit_z0 = rfittbl[:,4];
    rfit_p1p = rfittbl[:,5];
    rfit_p2p = rfittbl[:,6];
    rfit_p3p = rfittbl[:,7];
    rfit_p4p = rfittbl[:,8];
    rfit_chi2p = rfittbl[:,9];
    rfit_p1f = rfittbl[:,10];
    rfit_p2f = rfittbl[:,11];
    rfit_p3f = rfittbl[:,12];
    rfit_p4f = rfittbl[:,13];
    rfit_chi2f = rfittbl[:,14];
    rfit_cfxy = rfittbl[:,15];
    rfit_cftxy = rfittbl[:,16];

    # Read in the segment data for the forward fit.
    # segID nPts chi2avg chi2min chi2max
    segtbl = np.loadtxt("{0}/{1}/seg_{2}_{3}.dat".format(fit_outdir,trk_name,trk_name,ntrk));
    
    # Ensure there are multiple lines.
    if(len(segtbl.shape) == 1):
        seg_ID = []; seg_ID.append(segtbl[0]);
        seg_npts = []; seg_npts.append(segtbl[1]);
        seg_chi2avg = []; seg_chi2avg.append(segtbl[2]);
        seg_chi2min = []; seg_chi2min.append(segtbl[3]);
        seg_chi2max = []; seg_chi2max.append(segtbl[4]);
    else:
        seg_ID = segtbl[:,0];
        seg_npts = segtbl[:,1];
        seg_chi2avg = segtbl[:,2];
        seg_chi2min = segtbl[:,3];
        seg_chi2max = segtbl[:,4];
    
    # Read in the segment data for the reverse fit.
    # segID nPts chi2avg chi2min chi2max
    rsegtbl = np.loadtxt("{0}/{1}/rev/seg_{2}_{3}.dat".format(fit_outdir,trk_name,trk_name,ntrk));
    if(len(rsegtbl.shape) == 1):
        rseg_ID = []; rseg_ID.append(rsegtbl[0]);
        rseg_npts = []; rseg_npts.append(rsegtbl[1]);
        rseg_chi2avg = []; rseg_chi2avg.append(rsegtbl[2]);
        rseg_chi2min = []; rseg_chi2min.append(rsegtbl[3]);
        rseg_chi2max = []; rseg_chi2max.append(rsegtbl[4]);
    else:
        rseg_ID = rsegtbl[:,0];
        rseg_npts = rsegtbl[:,1];
        rseg_chi2avg = rsegtbl[:,2];
        rseg_chi2min = rsegtbl[:,3];
        rseg_chi2max = rsegtbl[:,4];
    
    if(len(seg_ID) == 0 or len(rseg_ID) == 0):
        print "ERROR track {0}: no segments exist in forward and/or reverse track".format(ntrk);
        continue;
    elif(len(seg_ID) < stat_nseg or len(rseg_ID) < stat_nseg):
        print "WARNING track {0}: fewer than stat_nseg = {1} segments exist in forward and/or reverse track".format(ntrk,stat_nseg);
    
    # Declare the relevant quantities.
    nkpts = 0; rnkpts = 0;
    akpts = 0; rakpts = 0;  # k for the purposes of calculating the average
    achi2 = 0.; rachi2 = 0.;
    
    # Extract the relevant quantities from the forward track.
    ss = 0;
    while(ss < stat_nseg and ss < len(seg_ID)):
        if(seg_chi2avg[ss] > 0):
            akpts += seg_npts[ss];
            achi2 += seg_npts[ss]*seg_chi2avg[ss];
        nkpts += seg_npts[ss];
        ss += 1;
    if(akpts > 0):
        achi2 /= akpts;
    else:
        achi2 = -1.;
        #print "WARNING: track {0} had a poor forward fit.".format(ntrk);
    
    splot_ktot.append(nkpts);
    splot_chi2avg.append(achi2);

    # Find the number of segments in the first fraction of the points (forward track).    
    ss = 0;
    nfrac = -1; nkpts = 0;
    ntotpts = sum(seg_npts);
    while(nfrac < 0 and ss < len(seg_ID)):
        nkpts += seg_npts[ss];
        if(nkpts > ntotpts*stat_efac):
            nfrac = seg_ID[ss]+1;
        ss += 1;

    nseg_ff = nfrac;
    nseg_fl = (seg_ID[-1]+1-nfrac);
    
    splot_sfac.append(1.0*nseg_fl/nseg_ff);
    
    # Compute the number of points per segment (forward track).
    if(nseg_ff > 0): splot_ptps_ff.append(int(ntotpts*stat_efac)/nseg_ff);
    else: splot_ptps_ff.append(-1.);
    if(nseg_fl > 0): splot_ptps_fl.append((ntotpts-int(ntotpts*stat_efac))/nseg_fl);
    else: splot_ptps_fl.append(-1.);
    
    # Compute the average chi2 for the first fraction of the points (forward track).
    nn = 0;
    nkapts = 0; chi2avg = 0.;
    nkacfxypts = 0; cfxyavg = 0.;
    while(nn < ntotpts*stat_efac):
        if(fit_chi2f[nn] > chi2_low and fit_chi2f[nn] < chi2_outlier):
            nkapts += 1;
            chi2avg += fit_chi2f[nn];
        if(fit_cfxy[nn] > cfxy_low and fit_cfxy[nn] < cfxy_outlier):
            nkacfxypts += 1;
            cfxyavg += fit_cfxy[nn];
        nn += 1;
    if(nkapts > 0): splot_schi2_ff.append(chi2avg/nkapts);
    else: splot_schi2_ff.append(-1.);
    if(nkacfxypts > 0): splot_scfxy_ff.append(cfxyavg/nkacfxypts);
    else: splot_scfxy_ff.append(-1.);
    #print "FF: Added avg. cfxy = {0}".format(cfxyavg/nkacfxypts)
    
    # Compute the average chi2 for the last fraction of the points (forward track).
    nkapts = 0; chi2avg = 0;
    nkacfxypts = 0; cfxyavg = 0.;
    while(nn < ntotpts):
        if(fit_chi2f[nn] > chi2_low and fit_chi2f[nn] < chi2_outlier):
            nkapts += 1;
            chi2avg += fit_chi2f[nn];
        if(fit_cfxy[nn] > cfxy_low and fit_cfxy[nn] < cfxy_outlier):
            nkacfxypts += 1;
            cfxyavg += fit_cfxy[nn];
        nn += 1;
    if(nkapts > 0): splot_schi2_fl.append(chi2avg/nkapts);
    else: splot_schi2_fl.append(-1.);
    if(nkacfxypts > 0): splot_scfxy_fl.append(cfxyavg/nkacfxypts);
    else: splot_scfxy_fl.append(-1.);
    #print "FL: Added avg. cfxy = {0}".format(cfxyavg/nkacfxypts)
        
    # Record the k/N and chi2 and cfxy values for points in the forward track
    #  with valid chi2.
    nn = 0;
    fchi2avg = 0.; navg = 0;
    cfxy_min = -1.; k_min = -1;
    while(nn < ntotpts):
        chi2 = fit_chi2f[nn];
        if(chi2 < chi2_outlier):
            prof_fkon.append(1.0*nn/ntotpts);
            prof_fchi2.append(chi2);
            fchi2avg += chi2; navg += 1
        cfxy = fit_cfxy[nn];
        if(cfxy < cfxy_outlier):
            prof_fckon.append(1.0*nn/ntotpts); 
            prof_fcfxy.append(cfxy);
            if(1.0*nn/ntotpts > 0.025 and (cfxy_min < 0 or cfxy < cfxy_min)):
                k_min = nn;
                cfxy_min = cfxy;
        nn += 1;
    if(fchi2avg/navg < 0.5): print "forward fchi2avg = {0}".format(fchi2avg/navg);
    splot_cfxy_min.append(cfxy_min);
    splot_cfxy_konmin.append(1.0*k_min/ntotpts);
    if(k_min < 0): print "\n\n** WARNING: kon_min < 0 for track {0}; ntotpts = {1}, k_min = {2}, cfxy_min = {3}".format(ntrk,ntotpts,k_min,cfxy_min);
 
    # Extract the relevant quantities from the reverse track.        
    ss = 0;
    while(ss < stat_nseg and ss < len(rseg_ID)):
        if(rseg_chi2avg[ss] > 0):
            rakpts += rseg_npts[ss];
            rachi2 += rseg_npts[ss]*rseg_chi2avg[ss];
        rnkpts += rseg_npts[ss];
        ss += 1;
    if(rakpts > 0):
        rachi2 /= rakpts;
    else:
        rachi2 = -1.;
        #print "WARNING: track {0} had a poor reverse fit.".format(ntrk);

    splot_rktot.append(rnkpts);
    splot_rchi2avg.append(rachi2);
    
    # Find the number of segments in the first fraction of the points (reverse track).   
    ss = 0;
    rnfrac = -1; rnkpts = 0;
    rntotpts = sum(rseg_npts);
    while(rnfrac < 0 and ss < len(rseg_ID)):
        rnkpts += rseg_npts[ss];
        if(rnkpts > rntotpts*stat_efac):
            rnfrac = rseg_ID[ss]+1;
        ss += 1;
    
    rnseg_ff = rnfrac;
    rnseg_fl = (rseg_ID[-1]+1-rnfrac);
    
    splot_rsfac.append(1.0*rnseg_fl/rnseg_ff);
    
    # Compute the number of points per segment (reverse track).
    if(rnseg_ff > 0): splot_ptps_rf.append(int(rntotpts*stat_efac)/rnseg_ff);
    else: splot_ptps_rf.append(-1.);
    if(rnseg_fl > 0): splot_ptps_rl.append((rntotpts-int(rntotpts*stat_efac))/rnseg_fl);
    else: splot_ptps_rl.append(-1.);

    # Compute the average chi2 and cfxy for the first fraction of the points (reverse track).
    nn = 0;
    rnkapts = 0; rchi2avg = 0.;
    rnkacfxypts = 0; rcfxyavg = 0.;
    while(nn < rntotpts*stat_efac):
        if(rfit_chi2f[nn] > chi2_low and rfit_chi2f[nn] < chi2_outlier):
            rnkapts += 1;
            rchi2avg += rfit_chi2f[nn];
        if(rfit_cfxy[nn] > cfxy_low and rfit_cfxy[nn] < cfxy_outlier):
            rnkacfxypts += 1;
            rcfxyavg += rfit_cfxy[nn];
        nn += 1;
    if(rnkapts > 0): splot_schi2_rf.append(rchi2avg/rnkapts);
    else: splot_schi2_rf.append(-1.);
    if(rnkacfxypts > 0): splot_scfxy_rf.append(rcfxyavg/rnkacfxypts);
    else: splot_scfxy_rf.append(-1.);
    #print "RF: Added avg. cfxy = {0}".format(rcfxyavg/rnkacfxypts)

    # Compute the average chi2 and cfxy for the last fraction of the points (reverse track).
    rnkapts = 0; rchi2avg = 0;
    rnkacfxypts = 0; rcfxyavg = 0.;
    while(nn < rntotpts):
        if(rfit_chi2f[nn] > chi2_low and rfit_chi2f[nn] < chi2_outlier):
            rnkapts += 1;
            rchi2avg += rfit_chi2f[nn];
        if(rfit_cfxy[nn] > cfxy_low and rfit_cfxy[nn] < cfxy_outlier):
            rnkacfxypts += 1;
            rcfxyavg += rfit_cfxy[nn];
        nn += 1;
    if(rnkapts > 0): splot_schi2_rl.append(rchi2avg/rnkapts);
    else: splot_schi2_rl.append(-1.);
    if(rnkacfxypts > 0): splot_scfxy_rl.append(rcfxyavg/rnkacfxypts);
    else: splot_scfxy_rl.append(-1.);
    #print "RL: Added avg. cfxy = {0}".format(rcfxyavg/rnkacfxypts)
        
    # Record the k/N and chi2 average for points in the reverse track
    #  with valid chi2.
    rnn = 0;
    rchi2avg = 0.; rnavg = 0;
    rcfxy_min = -1.; rk_min = -1;
    while(rnn < rntotpts):
        rchi2 = rfit_chi2f[rnn];
        if(rchi2 < chi2_outlier):
            prof_rkon.append(1.0*rnn/rntotpts);
            prof_rchi2.append(rchi2);
            rchi2avg += rchi2; rnavg += 1
        rcfxy = rfit_cfxy[rnn];
        if(rcfxy < cfxy_outlier):
            prof_rckon.append(1.0*rnn/rntotpts);
            prof_rcfxy.append(rcfxy);
            if(1.0*rnn/rntotpts > 0.025 and (rcfxy_min < 0 or rcfxy < rcfxy_min)):
                rk_min = rnn;
                rcfxy_min = rcfxy;
        rnn += 1;
    if(rchi2avg/rnavg < 0.5): print "reverse chi2 average is {0}".format(np.mean(rchi2avg/rnavg))
    splot_rcfxy_min.append(rcfxy_min);
    splot_rcfxy_konmin.append(1.0*rk_min/rntotpts);
    if(rk_min < 0): print "\n\n** WARNING: kon_min < 0 for track {0}".format(ntrk);

# ---------------------------------------------------------------------------
# Create arrays with the (last)/(first) chi2 averages.
splot_lof_forward = [];
for ff,ll in zip(splot_schi2_ff,splot_schi2_fl):
    splot_lof_forward.append(ll/ff);
splot_lof_reverse = [];
for ff,ll in zip(splot_schi2_rf,splot_schi2_rl):
    splot_lof_reverse.append(ll/ff);
    
# Create arrays with the (last)/(first) number of points / segment.
splot_lofptps_forward = [];
for ff,ll in zip(splot_ptps_ff,splot_ptps_fl):
    splot_lofptps_forward.append(ll/ff);
splot_lofptps_reverse = [];
for ff,ll in zip(splot_ptps_rf,splot_ptps_rl):
    splot_lofptps_reverse.append(ll/ff);

# ---------------------------------------------------------------------------
# Place the cut in k vs. chi2 for stat_nseg segments.
npass = 0; nrpass = 0;
for k,chi2,rk,rchi2 in zip(splot_ktot,splot_chi2avg,splot_rktot,splot_rchi2avg):
    
    # Check if the forward fit event has passed the cut.
    if(chi2 > 0 and chi2 < mcut*k + bcut):
        npass += 1;
    
    # Check if the reverse fit event has passed the cut.
    if(rchi2 > 0 and rchi2 < mcut*rk + bcut):
        nrpass += 1;
print "k vs. chi2 cut passes {0}% of forward fits and {1}% of reverse fits".format(num_tracks*npass/100.,num_tracks*nrpass/100.);

# Place the cut in s vs. chi2_avg.
npass = 0; nrpass = 0;
for s,loff,rs,lofr in zip(splot_sfac,splot_lof_forward,splot_rsfac,splot_lof_reverse):
    
    # Check if the forward fit event has passed the cut.
    if(loff > 0 and loff > smcut*s + sbcut):
        npass += 1;
    
    # Check if the reverse fit event has passed the cut.
    if(lofr > 0 and lofr > smcut*rs + sbcut):
        nrpass += 1;
print "s vs. chi2_avg cut passes {0}% of forward fits and {1}% of reverse fits".format(100.*npass/num_tracks,100.*nrpass/num_tracks);

# ---------------------------------------------------------------------------
# Determine the efficiency vs. 1-b using cfxy ratios.
fratio_cfxy = []; rratio_cfxy = [];
diffratio_cfxy = [];
for ffval, flval, rfval, rlval in zip(splot_scfxy_ff,splot_scfxy_fl,splot_scfxy_rf,splot_scfxy_rl):
    if(flval > 0. and ffval > 0. and rlval > 0. and rfval > 0.): 
        fratio_cfxy.append(ffval/flval);
        rratio_cfxy.append(rfval/rlval);
        diffratio_cfxy.append(rfval/rlval - ffval/flval);

eff_vals = []; bgr_vals = [];
cut_vals = [];
c1max = max(fratio_cfxy);
c2max = max(rratio_cfxy);
cmax = max(c1max,c2max);
nrvals = len(splot_scfxy_ff);
for cval in np.arange(0,cmax,cmax/100.):

    # Determine the number of signal and background events passing the cuts.
    nsig = 0; nbg = 0;
    for fr,rr in zip(fratio_cfxy,rratio_cfxy):
        if(fr < cval): nsig += 1;
        if(rr < cval): nbg += 1;

    # Add the fractions of signal and background events to the list.
    eff_vals.append(1.0*nsig/nrvals);
    bgr_vals.append(1.0-1.0*nbg/nrvals);

    # Add the cut to the cut values list.
    cut_vals.append(cval);

# ---------------------------------------------------------------------------
# Create the chi2 and cfxy vs. k/N profiles.
prof_kon = [];
fprof_nvals = []; fprof_chi2 = []; fprof_sigma = [];
rprof_nvals = []; rprof_chi2 = []; rprof_sigma = [];

fcprof_nvals = []; fcprof_cfxy = []; fcprof_sigma = [];
rcprof_nvals = []; rcprof_cfxy = []; rcprof_sigma = [];
for nn in range(nbins_kon):
    prof_kon.append(1.0*nn/nbins_kon);

    fprof_nvals.append(0); fprof_chi2.append(0.); fprof_sigma.append(0.);
    rprof_nvals.append(0); rprof_chi2.append(0.); rprof_sigma.append(0.);

    fcprof_nvals.append(0); fcprof_cfxy.append(0.); fcprof_sigma.append(0.);
    rcprof_nvals.append(0); rcprof_cfxy.append(0.); rcprof_sigma.append(0.);

for kon,fchi2 in zip(prof_fkon,prof_fchi2):
    bb = int(kon*nbins_kon);
    fprof_nvals[bb] += 1;
    fprof_chi2[bb] += fchi2;
    fprof_sigma[bb] += fchi2**2;

for kon,rchi2 in zip(prof_rkon,prof_rchi2):
    bb = int(kon*nbins_kon);
    rprof_nvals[bb] += 1;
    rprof_chi2[bb] += rchi2;
    rprof_sigma[bb] += rchi2**2;

for kon,fcfxy in zip(prof_fckon,prof_fcfxy):
    bb = int(kon*nbins_kon);
    fcprof_nvals[bb] += 1;
    fcprof_cfxy[bb] += fcfxy;
    fcprof_sigma[bb] += fcfxy**2;

for kon,rcfxy in zip(prof_rckon,prof_rcfxy):
    bb = int(kon*nbins_kon);
    rcprof_nvals[bb] += 1;
    rcprof_cfxy[bb] += rcfxy;
    rcprof_sigma[bb] += rcfxy**2;

# Normalize.
for bb in range(nbins_kon):
    if(fprof_nvals[bb] > 1):
        NN = fprof_nvals[bb];
        mu = fprof_chi2[bb]/fprof_nvals[bb];
        fprof_chi2[bb] = mu;
        fprof_sigma[bb] = sqrt((fprof_sigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);
    if(rprof_nvals[bb] > 1):
        NN = rprof_nvals[bb];
        mu = rprof_chi2[bb]/rprof_nvals[bb];
        rprof_chi2[bb] = mu;
        rprof_sigma[bb] = sqrt((rprof_sigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);
    if(fcprof_nvals[bb] > 1):
        NN = fcprof_nvals[bb];
        mu = fcprof_cfxy[bb]/fcprof_nvals[bb];
        fcprof_cfxy[bb] = mu;
        fcprof_sigma[bb] = sqrt((fcprof_sigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);
    if(rcprof_nvals[bb] > 1):
        NN = rcprof_nvals[bb];
        mu = rcprof_cfxy[bb]/rcprof_nvals[bb];
        rcprof_cfxy[bb] = mu;
        rcprof_sigma[bb] = sqrt((rcprof_sigma[bb])/(NN-1) - NN*mu**2/(NN-1))/sqrt(NN);

# Determine the number of events passed by different cuts in the cfxy minima.
cfxy_cuts = [];
cfxy_rnevts = []; cfxy_enevts = [];
cmax = 0.5;
for cval in np.arange(0,cmax,cmax/100.):
    
    # Count the number of events that pass this cut.
    rnpass = 0; enpass = 0;
    for min1,min2 in zip(splot_cfxy_konmin,splot_rcfxy_konmin):
        # Efficiency-biased cut
        if((min1 < cval and min2 > (1-cval)) or (min2 < cval and min1 > (1-cval))):
            enpass += 1;
        # Rejection-biased cut
        if(min1 < cval or min2 < cval or min1 > (1-cval) or min2 > (1-cval)):
            rnpass += 1;
    
    cfxy_cuts.append(cval);
    cfxy_rnevts.append(rnpass); cfxy_enevts.append(enpass);

# Print out the profiles.
f_fprof = open("{0}/prof/prof_{1}_forward.dat".format(fit_outdir,trk_name),"w")
f_fprof.write("# (k/N) (chi2) (sigma) (N)\n")
for kon,chi2,sigma,NN in zip(prof_kon,fprof_chi2,fprof_sigma,fprof_nvals):
    f_fprof.write("{0} {1} {2} {3}\n".format(kon,chi2,sigma,NN));
f_fprof.close()

f_rprof = open("{0}/prof/prof_{1}_reverse.dat".format(fit_outdir,trk_name),"w")
f_rprof.write("# (k/N) (chi2) (sigma) (N)\n")
for kon,chi2,sigma,NN in zip(prof_kon,rprof_chi2,rprof_sigma,rprof_nvals):
    f_rprof.write("{0} {1} {2} {3}\n".format(kon,chi2,sigma,NN));
f_rprof.close()

f_fcprof = open("{0}/prof/cprof_{1}_forward.dat".format(fit_outdir,trk_name),"w")
f_fcprof.write("# (k/N) (cfxy) (sigma) (N)\n")
for kon,cfxy,sigma,NN in zip(prof_kon,fcprof_cfxy,fcprof_sigma,fcprof_nvals):
    f_fcprof.write("{0} {1} {2} {3}\n".format(kon,cfxy,sigma,NN));
f_fcprof.close()

f_rcprof = open("{0}/prof/cprof_{1}_reverse.dat".format(fit_outdir,trk_name),"w")
f_rcprof.write("# (k/N) (cfxy) (sigma) (N)\n")
for kon,cfxy,sigma,NN in zip(prof_kon,rcprof_cfxy,rcprof_sigma,rcprof_nvals):
    f_rcprof.write("{0} {1} {2} {3}\n".format(kon,cfxy,sigma,NN));
f_rcprof.close()

# Print the number of passing events vs. cut value.
f_cfxycuts = open("{0}/cfxy_cuts.dat".format(plt_base),"w")
f_cfxycuts.write("# (cut value) (pass fraction rej. biased) (pass fraction eff. biased)\n")
for cval,nrpass,nepass in zip(cfxy_cuts,cfxy_rnevts,cfxy_enevts):
    f_cfxycuts.write("{0} {1} {2}\n".format(cval,1.0*nrpass/num_tracks,1.0*nepass/num_tracks));
f_cfxycuts.close();

# ---------------------------------------------------------------------------
# Prepare to plot the k vs. chi2 cut.
xplot = np.linspace(min(splot_ktot),max(splot_ktot),100);
zcut = [mcut, bcut];
pcut = np.poly1d(zcut);

# Prepare to plot the s vs. chi2_avg cut.
xsplot = np.linspace(min(splot_sfac),max(splot_sfac),100);
zscut = [smcut, sbcut];
pscut = np.poly1d(zscut);

# Plot the quantities in a scatter plot.
fig = plt.figure(1);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_ktot,splot_chi2avg,'.',color='blue',markersize=1.0,label='Forward fit');
plt.plot(splot_rktot,splot_rchi2avg,'.',color='red',markersize=1.0,label='Reverse fit');
plt.plot(xplot,pcut(xplot),'--',color='black',linewidth=1);
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.title("First two segments ($\chi^{2}_{\mathrm{limit}}$ = 20); avg. $\chi^{2}$ vs. k$_{\mathrm{total}}$");
plt.xlabel("k$_{\mathrm{total}}$");
plt.ylabel("$\chi^{2}_{\mathrm{avg}}$");
plt.savefig("{0}/seg_kvschi2.pdf".format(plt_base), bbox_inches='tight');

# Plot the chi2-values in a histogram.
fig = plt.figure(2);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
chi2n, chi2bins, chi2patches = plt.hist(splot_chi2avg, 20, normed=0, histtype='step',color='blue',label='Forward fit');
rchi2n, rchi2bins, rchi2patches = plt.hist(splot_rchi2avg, 20, normed=0, histtype='step',color='red',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("$\chi^{2}_{\mathrm{avg}}$");
plt.ylabel("Counts/bin");
plt.savefig("{0}/seg_chi2hist.pdf".format(plt_base), bbox_inches='tight');

# Plot the k-values in a histogram.
fig = plt.figure(3);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(splot_ktot, 50, normed=0, histtype='step',color='blue',label='Forward fit');
krn, krbins, krpatches = plt.hist(splot_rktot, 50, normed=0, histtype='step',color='red',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.axis([0,max(splot_ktot),0,max(kn)]);
plt.xlabel("k$_{\mathrm{total}}$");
plt.ylabel("Counts/bin");
plt.savefig("{0}/seg_khist.pdf".format(plt_base), bbox_inches='tight');

# Plot the s-values in a histogram.
fig = plt.figure(4);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
sn, sbins, spatches = plt.hist(splot_sfac, 20, normed=0, histtype='step',color='blue',label='Forward fit');
srn, srbins, srpatches = plt.hist(splot_rsfac, 20, normed=0, histtype='step',color='red',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("s");
plt.ylabel("Counts/bin");
plt.savefig("{0}/seg_svalhist.pdf".format(plt_base), bbox_inches='tight');

# Scatter plot of the average chi2 for the forward and reverse tracks.
fig = plt.figure(5);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_schi2_ff, splot_schi2_fl, '.', color='red',label='Forward fit');
plt.plot(splot_schi2_rf, splot_schi2_rl, '.', color='blue',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("$\chi^{2}$ (first)");
plt.ylabel("$\chi^{2}$ (last)");
plt.savefig("{0}/seg_chi2avg.pdf".format(plt_base), bbox_inches='tight');

# Scatter plot of the average cfxy value for the forward and reverse tracks.
fig = plt.figure(6);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_scfxy_ff, splot_scfxy_fl, '.', color='red');
plt.plot(splot_scfxy_rf, splot_scfxy_rl, '.', color='blue');
plt.xlabel("C_F,xy (first)");
plt.ylabel("C_F,xy (last)");
plt.savefig("{0}/scatter_cfxyavg.pdf".format(plt_base), bbox_inches='tight');

# Histogram of the average cfxy ratios.
fig = plt.figure(7);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
sn, sbins, spatches = plt.hist(fratio_cfxy, 100, normed=0, histtype='step',color='blue',label='Forward fit');
srn, srbins, srpatches = plt.hist(rratio_cfxy, 100, normed=0, histtype='step',color='red',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("C_F,xy (first) / C_F,xy (last)");
plt.ylabel("Counts/bin");
plt.savefig("{0}/hist_cfxyavg.pdf".format(plt_base), bbox_inches='tight');

# Histogram the difference in the ratios. 
fig = plt.figure(8);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
sn, sbins, spatches = plt.hist(diffratio_cfxy, 100, normed=0, histtype='step',color='blue');
plt.xlabel("C_F,xy (first) / C_F,xy (last) difference (reverse - forward)");
plt.ylabel("Counts/bin");
plt.savefig("{0}/hist_rdiffcfxy.pdf".format(plt_base), bbox_inches='tight');

# Scatter plot of s-value vs. last/first chi2 averages.
fig = plt.figure(9);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_sfac, splot_lof_forward, '.', color='red', label='Forward fit');
plt.plot(splot_rsfac, splot_lof_reverse, '.', color='blue', label='Reverse fit');
plt.plot(xsplot,pscut(xsplot),'--',color='black',linewidth=1);
plt.title("s vs. chi2_avg ratio: f = {0}, $\chi^2$ limit = {1}".format(stat_efac,chi2_outlier));
lnd = plt.legend(loc=3,frameon=False,handletextpad=0);
plt.xlabel("s");
plt.ylabel("avg. ($\chi^{2}$ last/$\chi^{2}$ first)");
plt.savefig("{0}/seg_s_vs_chi2avg.pdf".format(plt_base), bbox_inches='tight');

# Scatter plot of s-value vs. number of points/segment ratios.
fig = plt.figure(10);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(splot_sfac, splot_lofptps_forward, '.', color='red', label='Forward fit');
plt.plot(splot_rsfac, splot_lofptps_reverse, '.', color='blue', label='Reverse fit');
plt.plot(xsplot,pscut(xsplot),'--',color='black',linewidth=1);
plt.title("s vs. points/segment ratio: f = {0}, $\chi^2$ limit = {1}".format(stat_efac,chi2_outlier));
lnd = plt.legend(loc=3,frameon=False,handletextpad=0);
plt.xlabel("s");
plt.ylabel("avg. ($\chi^{2}$ last/$\chi^{2}$ first)");
plt.savefig("{0}/seg_s_vs_ptps.pdf".format(plt_base), bbox_inches='tight');

# Plot and print the chi2 profiles.
fig = plt.figure(11);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, fprof_chi2, '.', color='red', label='Forward fit');
plt.plot(prof_kon, rprof_chi2, '.', color='blue', label='Reverse fit');
plt.title("chi2 profile: $\chi^2$ limit = {0}".format(chi2_outlier));
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("Track fraction (k/N)");
plt.ylabel("$\chi^{2}$/dof average");
plt.savefig("{0}/fit_profiles.pdf".format(plt_base), bbox_inches='tight');

# Plot and print the cfxy profiles.
fig = plt.figure(12);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(prof_kon, fcprof_cfxy, '.', color='red', label='Forward fit');
plt.plot(prof_kon, rcprof_cfxy, '.', color='blue', label='Reverse fit');
plt.title("cfxy profile: cfxy_limit = {0}".format(cfxy_outlier));
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("Track fraction (k/N)");
plt.ylabel("cfxy average");
plt.savefig("{0}/cfxy_profiles.pdf".format(plt_base), bbox_inches='tight');

# Histogram the cfxy minimum k/N values.
fig = plt.figure(13);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
kn, kbins, kpatches = plt.hist(splot_cfxy_konmin, 100, normed=0, histtype='step',color='blue',label='Forward fit');
krn, krbins, krpatches = plt.hist(splot_rcfxy_konmin, 100, normed=0, histtype='step',color='red',label='Reverse fit');
lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
plt.xlabel("C_F,xy k/N at minimum");
plt.ylabel("Counts/bin");
plt.savefig("{0}/hist_cfxy_konmin.pdf".format(plt_base), bbox_inches='tight');

# Plot the signal eff. vs. background rejection. 
fig = plt.figure(14);
fig.set_figheight(5.0);
fig.set_figwidth(7.5);
plt.plot(bgr_vals,eff_vals,'-',color='blue',markersize=1.0);
#lnd = plt.legend(loc=1,frameon=False,handletextpad=0);
#plt.axis([0,max(splot_ktot),0,max(splot_fratio)]);
plt.xlabel("Background rejection (1-b)");
plt.ylabel("Signal efficiency ($\epsilon$)");
plt.savefig("{0}/cfxy_signal_eff_vs_bgr.pdf".format(plt_base), bbox_inches='tight');
