"""
trackgen.py

Generates tracks with multiple scattering.

Track output format:
    
    x0 y0 z0 zi zf ux uy uz E deltaE deltaX
    
    (x0,y0): what would be the recorded slice (x,y);
    zi, zf: initial and final z locations
    ux,uy,uz: direction vector for the slcie
    E: energy of the track at this point
    deltaE: the energy deposited in this slice
    deltaX: the step size for this slice
    
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

from abc import ABCMeta, abstractmethod
import logging 
logging.basicConfig(level=logging.INFO)

if(not os.path.isdir(trk_outdir)): os.mkdir(trk_outdir);
if(not os.path.isdir("{0}/{1}")): os.mkdir(trk_outdir,trk_name);
if(not os.path.isdir("{0}/{1}/rev".format(trk_outdir,trk_name))): os.mkdir("{0}/{1}/rev".format(trk_outdir,trk_name));

# Read in the stopping power and interpolate.
xesp_tbl = np.loadtxt("data/xe_estopping_power_NIST.dat");
rho = pc_rho0*(Pgas/(Tgas/273.15))*(pc_m_Xe/pc_NA);
e_vals = xesp_tbl[:,0];
dEdx_vals = xesp_tbl[:,1];
e_vals = np.insert(e_vals,0,0.0);
dEdx_vals = np.insert(dEdx_vals,0,dEdx_vals[0]);
xesp = interp1d(e_vals,dEdx_vals*rho,kind='cubic');
print e_vals;
print dEdx_vals;

# Create num_tracks tracks.
for ntrk in range(num_tracks):
    
    # Declare arrays for the track.
    # x0 y0 zi zf ux uy uz E deltaE deltaX
    trk_x = []; trk_y = []; trk_zi = []; trk_zf = [];
    trk_ux = []; trk_uy = []; trk_uz = [];
    trk_E = []; trk_deltaE = []; trk_deltaX = [];

    # Initialize the track.
    te = E_0;
    tx = 0.; ty = 0.; tz = 0.;
    ux = 0.; uy = 0.; uz = 1.;
    
    logging.debug("\n\n-- Track {0} --\n\n".format(ntrk));
    
    # Continue until 0 energy.
    while(te > E_tol):

        # Compute the current momentum of the track.
        ptrk = sqrt((te + 0.511)**2 - 0.511**2); 
        
        # Determine the energy loss for this step.        
        if(te < eslice):
            deltaE = te;
        else:
            deltaE = eslice;
        te -= deltaE;
        if(te < 0.): te = 0.;
        
        logging.debug("-> Energy = {0}, energy loss: {1}".format(te,deltaE));

        # Determine the distance of this step.
        deltaX = integrate.quad(lambda ee: 1./xesp(ee),te,te+deltaE,limit=1000)[0]
        
        # Make the step.
        dx = deltaX*ux; dy = deltaX*uy; dz = deltaX*uz;
        
        logging.debug("-> Step: deltaX = {0}; dx = {1}, dy = {2}, dz = {3}".format(deltaX,dx,dy,dz));
        
        # Record the variables for the step.
        trk_x.append(tx + dx/2.);
        trk_y.append(ty + dy/2.);
        trk_zi.append(tz);
        trk_zf.append(tz + dz);
        trk_ux.append(ux);
        trk_uy.append(uy);
        trk_uz.append(uz);
        trk_E.append(te + deltaE);
        trk_deltaE.append(deltaE);
        trk_deltaX.append(deltaX);
        
        # Update the positions.
        tx += dx; ty += dy; tz += dz;
        
        # Determine the scattering angles in the frame in which the track
        #  direction is aligned with the z-axis.
        sigma_theta = SigmaThetaMs(ptrk,deltaX/Lr);
        if(deltaX/Lr < 0.001 or deltaX/Lr > 100.):
            print "WARNING: L/Lr = {0} out of range of validity of the formula.".format(deltaX/Lr);
        thetaX = rd.gauss(0,sigma_theta);
        thetaY = rd.gauss(0,sigma_theta);
        tanX = tan(thetaX);
        tanY = tan(thetaY);
        
        logging.debug("-> sigma(theta) = {0}; tanX = {1}, tanY = {2}".format(sigma_theta,tanX,tanY));
        
        # Compute the direction cosines of the rotation matrix to move to the lab frame.
        nxy = sqrt(ux**2 + uy**2);
        if(nxy > 0.):
            alpha1 = uy/nxy; alpha2 = -ux*uz/nxy; alpha3 = ux;
            beta1 = -ux/nxy; beta2 = -uy*uz/nxy; beta3 = uy;
            gamma1 = 0.; gamma2 = nxy; gamma3 = uz;
        else:
            # Special case; the direction vector is the x-axis.  Choose
            #  the orthonormal basis as the normal unit vectors.
            alpha1 = 1.; alpha2 = 0.; alpha3 = 0.;
            beta1 = 0.; beta2 = 1.; beta3 = 0.;
            gamma1 = 0.; gamma2 = 0.; gamma3 = 1.;
        
        logging.debug("-> alpha1 = {0}, alpha2 = {1}, alpha3 = {2}".format(alpha1,alpha2,alpha3));
        logging.debug("-> alpha1 = {0}, alpha2 = {1}, alpha3 = {2}".format(alpha1,alpha2,alpha3));
        logging.debug("-> beta1 = {0}, beta2 = {1}, beta3 = {2}".format(beta1,beta2,beta3));
        logging.debug("-> gamma1 = {0}, gamma2 = {1}, gamma3 = {2}".format(gamma1,gamma2,gamma3));

        # Determine direction vector components in the reference (lab) frame.
        nrm = sqrt(tanX**2 + tanY**2 + 1);
        xp = (alpha1*tanX + alpha2*tanY + alpha3)/nrm;
        yp = (beta1*tanX + beta2*tanY + beta3)/nrm;
        zp = (gamma1*tanX + gamma2*tanY + gamma3)/nrm;
        
        # Set the new direction vector.
        nrm = sqrt(xp**2 + yp**2 + zp**2);
        ux = xp/nrm;
        uy = yp/nrm;
        uz = zp/nrm;
    
    # Print out the track.
    fn_trk = "{0}/{1}/{2}_{3}.dat".format(trk_outdir,trk_name,trk_name,ntrk);
    f_trk = open(fn_trk,"w");
    f_trk.write("# x0 y0 zi zf ux uy uz E deltaE deltaX\n");
    for x0_f,y0_f,zi_f,zf_f,ux_f,uy_f,uz_f,E_f,deltaE_f,deltaX_f in zip(trk_x,trk_y,trk_zi,trk_zf,trk_ux,trk_uy,trk_uz,trk_E,trk_deltaE,trk_deltaX):
        f_trk.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(x0_f,y0_f,zi_f,zf_f,ux_f,uy_f,uz_f,E_f,deltaE_f,deltaX_f));
    f_trk.close();
    
    # -----------------------------------------------------------------------
    # Create the reversed track.
    # -----------------------------------------------------------------------
    
    # Reverse all lists.
    trk_x.reverse(); trk_y.reverse(); trk_zi.reverse(); trk_zf.reverse();
    trk_ux.reverse(); trk_uy.reverse(); trk_uz.reverse();
    trk_E.reverse(); trk_deltaE.reverse(); trk_deltaX.reverse();
    
    # Get the initial direction.
    zi0 = (trk_zi[0] + trk_zf[0])/2.;
    zi1 = (trk_zi[1] + trk_zf[1])/2.;
    ux = (trk_x[1] - trk_x[0]);
    uy = (trk_y[1] - trk_y[0]);
    uz = (zi1 - zi0);
    imag = sqrt(ux**2 + uy**2 + uz**2);
    ux /= imag;
    uy /= imag;
    uz /= imag;
    
#    # Rotate the entire track such that the initial direction is the z-direction.
#    nxy = sqrt(ux**2 + uy**2);
#    if(nxy > 0.):
#        alpha1 = uy/nxy; alpha2 = -ux*uz/nxy; alpha3 = ux;
#        beta1 = -ux/nxy; beta2 = -uy*uz/nxy; beta3 = uy;
#        gamma1 = 0.; gamma2 = nxy; gamma3 = uz;
#    else:
#        # Special case; the direction vector is the x-axis.  Choose
#        #  the orthonormal basis as the normal unit vectors.
#        alpha1 = 1.; alpha2 = 0.; alpha3 = 0.;
#        beta1 = 0.; beta2 = 1.; beta3 = 0.;
#        gamma1 = 0.; gamma2 = 0.; gamma3 = 1.;
    
    # Shift the origin and rotate the track.
    rtrk_x = []; rtrk_y = []; rtrk_zi = []; rtrk_zf = [];
    rtrk_ux = []; rtrk_uy = []; rtrk_uz = [];
    rtrk_E = []; rtrk_deltaE = []; rtrk_deltaX = [];
    tx0 = trk_x[0]; ty0 = trk_y[0]; tz0 = trk_zf[0];
    for x0_f,y0_f,zi_f,zf_f,ux_f,uy_f,uz_f,E_f,deltaE_f,deltaX_f in zip(trk_x,trk_y,trk_zi,trk_zf,trk_ux,trk_uy,trk_uz,trk_E,trk_deltaE,trk_deltaX):
        
        # Shift the origin and reflect.
        x0 = -1*(x0_f - tx0);
        y0 = -1*(y0_f - ty0);
        zi0 = -1*(zi_f - tz0);
        zf0 = -1*(zf_f - tz0);
        z0 = (zi0 + zf0)/2.;
        
#        # Rotate.
#        xf = alpha1*x0 + alpha2*y0 + alpha3*z0;
#        yf = beta1*x0 + beta2*y0 + beta3*z0;
#        zif = gamma1*x0 + gamma2*y0 + gamma3*zi0;
#        zff = gamma1*x0 + gamma2*y0 + gamma3*zf0;
        
        # Store the shifted and reflected values: note flip zi <-> zf
        rtrk_x.append(x0);
        rtrk_y.append(y0);
        rtrk_zi.append(zf0);
        rtrk_zf.append(zi0);
        
        # For now keep the old direction vectors (this is not something we will have in the real data).
        rtrk_ux.append(ux_f);
        rtrk_uy.append(uy_f);
        rtrk_uz.append(uz_f);
        
        # Keep the old E, deltaE and deltaX.
        rtrk_E.append(E_f);
        rtrk_deltaE.append(deltaE_f);
        rtrk_deltaX.append(deltaX_f);

    # Print the reversed track.
    fn_rtrk = "{0}/{1}/rev/{2}_{3}.dat".format(trk_outdir,trk_name,trk_name,ntrk);
    f_rtrk = open(fn_rtrk,"w");
    f_rtrk.write("# x0 y0 zi zf ux uy uz E deltaE deltaX\n");
    for x0_f,y0_f,zi_f,zf_f,ux_f,uy_f,uz_f,E_f,deltaE_f,deltaX_f in zip(rtrk_x,rtrk_y,rtrk_zi,rtrk_zf,rtrk_ux,rtrk_uy,rtrk_uz,rtrk_E,trk_deltaE,trk_deltaX):
        f_rtrk.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(x0_f,y0_f,zi_f,zf_f,ux_f,uy_f,uz_f,E_f,deltaE_f,deltaX_f));
    f_rtrk.close();
