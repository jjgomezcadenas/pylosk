"""
trackplot.py

Plots tracks and their Kalman filter fit.

"""
import sys
import numpy as np
import scipy.integrate as integrate
import random as rd
import matplotlib.pyplot as plt
import os
from mpl_toolkits.mplot3d import Axes3D
from math import *
from scipy.interpolate import interp1d
from trackdefs import *

from abc import ABCMeta, abstractmethod
import logging 
logging.basicConfig(level=logging.DEBUG)

if(rev_trk):
    print "\n\n-- WORKING ON REVERSED TRACKS --\n\n";
    #fnb_trk = "{0}/rev/".format(trk_outdir);
    fnb_trk = "{0}/{1}".format(trk_outdir,trk_name);
    fnb_fit = "{0}/{1}/rev".format(fit_outdir,trk_name);
    fnb_plt = "{0}/{1}/rev".format(plt_outdir,trk_name);
else: 
    fnb_trk = "{0}/{1}".format(trk_outdir,trk_name);
    fnb_fit = "{0}/{1}".format(fit_outdir,trk_name);
    fnb_plt = "{0}/{1}".format(plt_outdir,trk_name);
    
if(not os.path.isdir(fnb_trk) and not plt_smearedHits): 
    print "ERROR: attempting to plot original hits and track directory {0} not available".format(fnb_trk);
    sys.exit();
if(not os.path.isdir(fnb_fit)): 
    print "ERROR: fit directory {0} not available".format(fnb_fit);
    sys.exit();

if(not os.path.isdir("{0}/{1}".format(plt_outdir,trk_name))):
    os.mkdir("{0}/{1}".format(plt_outdir,trk_name));
    print "Creating plot directory {0}/{1}...".format(plt_outdir,trk_name);
if(not os.path.isdir("{0}/{1}/rev".format(plt_outdir,trk_name))):
    os.mkdir("{0}/{1}/rev".format(plt_outdir,trk_name));
    print "Creating plot directory {0}/{1}/rev...".format(plt_outdir,trk_name);

# Keep a running list of the values of chi2.
chi2_list = [];

# Create num_tracks tracks.
for ntrk in range(num_tracks):
    
    logging.debug("-- Printing track {0}\n".format(ntrk));

    # Read in the track.
    # x0 y0 zi zf ux uy uz E deltaE deltaX
    if(not plt_smearedHits):
        trktbl = np.loadtxt("{0}/{1}_{2}.dat".format(fnb_trk,trk_name,ntrk));
        trk_x0 = trktbl[:,0];
        trk_y0 = trktbl[:,1];
        trk_zi = trktbl[:,2];
        trk_zf = trktbl[:,3];
        trk_ux = trktbl[:,4];
        trk_uy = trktbl[:,5];
        trk_uz = trktbl[:,6];
        trk_E = trktbl[:,7];
        trk_deltaE = trktbl[:,8];
        trk_deltaX = trktbl[:,9];
        
        trk_z0 = [];
        for zi,zf in zip(trk_zi,trk_zf): trk_z0.append((zi+zf)/2.);
    
    # Read in the fit.
    # segID k p1p p2p p3p p4p chi2p p1f p2f p3f p4f chi2f
    fittbl = np.loadtxt("{0}/fit_{1}_{2}.dat".format(fnb_fit,trk_name,ntrk));
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
 
    fit_px0 = []; fit_py0 = []; fit_pz0 = [];
    for z0,p1p,p2p,p3p,p4p in zip(fit_z0,fit_p1p,fit_p2p,fit_p3p,fit_p4p):
        x0 = p1p + z0*p3p;
        y0 = p2p + z0*p4p;
        
        fit_px0.append(x0);
        fit_py0.append(y0);
        fit_pz0.append(z0);
    
    fit_fx0 = []; fit_fy0 = []; fit_fz0 = [];
    for z0,p1f,p2f,p3f,p4f,chi2 in zip(fit_z0,fit_p1f,fit_p2f,fit_p3f,fit_p4f,fit_chi2f):
        x0 = p1f + z0*p3f;
        y0 = p2f + z0*p4f;
        
        fit_fx0.append(x0);
        fit_fy0.append(y0);
        fit_fz0.append(z0);
        
        # Add the chi2 if it is not an outlier.
        if(chi2 < chi2_outlier):
            chi2_list.append(chi2);

    
    # Plot the 3-D track plot with projections and chi2.
    if(plt_tracks):
        
        # Make the plot.
        fig = plt.figure();
        fig.set_figheight(15.0);
        fig.set_figwidth(10.0);
        
        # Create the 3D track plot.
        ax1 = fig.add_subplot(321, projection='3d');
        if(plt_smearedHits):
            ax1.plot(fit_x0,fit_y0,fit_z0,'o');
        else:
            ax1.plot(trk_x0,trk_y0,trk_z0,'o');
        if(plt_filtered): ax1.plot(fit_fx0,fit_fy0,fit_fz0,'.',color='black');
        if(plt_prediction): ax1.plot(fit_px0[1:],fit_py0[1:],fit_pz0[1:],'.',color='blue');
        ax1.set_xlabel("x ({0})".format(plt_units));
        ax1.set_ylabel("y ({0})".format(plt_units));
        ax1.set_zlabel("z ({0})".format(plt_units));
#        print fit_px0
#        print fit_py0
#        print fit_pz0
        
        # Create the x-y projection.
        ax2 = fig.add_subplot(322);
        if(plt_smearedHits):
            ax2.plot(fit_x0,fit_y0,'o');
        else:
            ax2.plot(trk_x0,trk_y0,'o');
        if(plt_filtered): ax2.plot(fit_fx0,fit_fy0,'.',color='black');
        if(plt_prediction): ax2.plot(fit_px0[1:],fit_py0[1:],'.',color='blue');
        ax2.set_xlabel("x ({0})".format(plt_units));
        ax2.set_ylabel("y ({0})".format(plt_units));
        
        # Create the x-z projection.
        ax3 = fig.add_subplot(323);
        if(plt_smearedHits):
            ax3.plot(fit_x0,fit_z0,'o');
        else:
            ax3.plot(trk_x0,trk_z0,'o');
        if(plt_filtered): ax3.plot(fit_fx0,fit_fz0,'.',color='black');
        if(plt_prediction): ax3.plot(fit_px0[1:],fit_pz0[1:],'.',color='blue');
        ax3.set_xlabel("x ({0})".format(plt_units));
        ax3.set_ylabel("z ({0})".format(plt_units));    
        
        # Create the y-z projection.
        ax4 = fig.add_subplot(324);
        if(plt_smearedHits):
            ax4.plot(fit_y0,fit_z0,'o');
        else:
            ax4.plot(trk_y0,trk_z0,'o');
        if(plt_filtered): ax4.plot(fit_fy0,fit_fz0,'.',color='black');
        if(plt_prediction): ax4.plot(fit_py0[1:],fit_pz0[1:],'.',color='blue');
        ax4.set_xlabel("y ({0})".format(plt_units));
        ax4.set_ylabel("z ({0})".format(plt_units));
        
        # Plot the chi2.
        ax5 = fig.add_subplot(325);
        ax5.plot(fit_k[1:],fit_chi2p[1:],color='black');
        ax5.plot(fit_k[1:],fit_chi2f[1:],color='blue');
        ax5.set_xlabel("k");
        ax5.set_ylabel("$\chi^{2}$");
        #ax5.set_yscale("log");
        ax5.set_title("Avg. chi2f = {0}".format(np.mean(fit_chi2f[1:])));

        # Plot the square root of the sum of the first two diagonals of the covariance matrix.
        cfxy_k = [];
        cfxy_vals = [];
        for k,cfxy in zip(fit_k[1:],fit_cfxy[1:]):
            if(cfxy > cfxy_low and cfxy < cfxy_outlier):
              cfxy_k.append(k);
              cfxy_vals.append(cfxy);
        ax6 = fig.add_subplot(326);
        ax6.plot(cfxy_k,cfxy_vals,color='black');
        #ax6.plot(fit_k[1:-1],np.diff(fit_cfxy[1:]),color='blue');
        ax6.set_xlabel("k");
        ax6.set_ylabel("C$_{F,xy}$");
        #ax6.set_yscale("log");
        #ax6.set_title("");       
 
        # Show and/or print the plot.
        if(plt_print):
            fn_plt = "{0}/plt_{1}_{2}.pdf".format(fnb_plt,trk_name,ntrk);
            plt.savefig(fn_plt, bbox_inches='tight');
        if(plt_show):
            plt.show();
            
# Plot the chi2 histogram.
if(plt_chi2):
    fig = plt.figure();
    fig.set_figheight(5.0);
    fig.set_figwidth(7.5);
    chi2n, chi2bins, chi2patches = plt.hist(chi2_list, 100, normed=0, histtype='step',color='blue',label='Forward fit');
    plt.xlabel("$\chi^{2}$");
    plt.ylabel("Counts/bin");
    
    # Show and/or print the plot.
    if(plt_print):
        fn_plt = "{0}/chi2_{1}.pdf".format(fnb_plt,trk_name);
        plt.savefig(fn_plt, bbox_inches='tight');
    if(plt_show):
        plt.show();
