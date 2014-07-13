"""
KTrackSegment.py

Contains state and chi2 information on a segment of the fitted track.

"""
import sys
import numpy as np
from KFBase import KFVector, KFMatrix
from KFMeasurement import KFMeasurement 
from KFZSlice import KFZSlice
from KFNode import KFNode 
from KFState import KFState
from KFSetup import KFSetup 
from KFSystem import KFSystem
from KalmanFilter import KalmanFilter
from math import *


class KTrackSegment(object):
    """
    A segment of a fitted track
    """
        
    def __init__(self,segID):
        """
        Initialize the track segment.
        """
         
        self.seg_id = segID;
        
        self.seg_k = []; self.seg_z0 = [];
        self.seg_p1p = []; self.seg_p2p = []; self.seg_p3p = []; self.seg_p4p = [];
        self.seg_p1f = []; self.seg_p2f = []; self.seg_p3f = []; self.seg_p4f = [];
        self.seg_cfxy = []; self.seg_cftxy = [];
        self.seg_pchisq = []; self.seg_fchisq = [];
        self.seg_x0 = []; self.seg_y0 = [];
        
    def AddPoint(self,k,x0,y0,z0,p1p,p2p,p3p,p4p,p1f,p2f,p3f,p4f,cfxy,cftxy,pchisq,fchisq):
        """
        Add a point to the segment.
        """
        
        self.seg_k.append(k);
        self.seg_z0.append(z0);
        self.seg_p1p.append(p1p); self.seg_p2p.append(p2p); 
        self.seg_p3p.append(p3p); self.seg_p4p.append(p4p);
        self.seg_p1f.append(p1f); self.seg_p2f.append(p2f); 
        self.seg_p3f.append(p3f); self.seg_p4f.append(p4f);
        self.seg_cfxy.append(cfxy); self.seg_cftxy.append(cftxy);
        self.seg_pchisq.append(pchisq); self.seg_fchisq.append(fchisq);
        self.seg_x0.append(x0); self.seg_y0.append(y0);

    def NumPoints(self):
        """
        Returns the number of points in the segment.
        """
        
        return len(seg_k)

    def PredictedState(self,i):
        """
        Returns the predicted state in point i of segment
        """
        PS=[]
        PS.append(self.seg_p1p[i])
        PS.append(self.seg_p2p[i])
        PS.append(self.seg_p3p[i])
        PS.append(self.seg_p4p[i])

        return PS

    def FilteredState(self,i):
        """
        Returns the filteres state in point i of segment
        """
        PS=[]
        PS.append(self.seg_p1f[i])
        PS.append(self.seg_p2f[i])
        PS.append(self.seg_p3f[i])
        PS.append(self.seg_p4f[i])

        return PS

    def PredictedChi2(self,i):
        """
        Returns the chi2 (predicted)
        """
        return self.seg_pchisq

    def FilteredChi2(self,i):
        """
        Returns the chi2 (filtered)
        """
        return self.seg_fchisq

        
