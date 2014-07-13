"""
A collection of utility functions
"""
from KLog import *
import KParam as kp
import numpy as np
import sys,os,getopt
from KalmanFilter import *
from IEventReader import IEventReader
from IParticle import IParticle
from TDrawHits import TDrawHits
from MPLDrawHits import MPLDrawHits
from MPLDrawVector import MPLDrawVector

debug = Debug.draw.value
DrawRoot = True
DrawMPL = True


def KSystem(Hits):
    """
    Returns a KFSystem  
    """    
    Nodes=[]    
    numberOfHits = len(Hits)
    Dz =[]
    Zi =[]
    for k in range(0,numberOfHits-1):
        hitk =Hits[k]
        hitkp1 =Hits[k+1]
        zki =hitk[2]
        zkf = hitkp1[2]
        edepk = hitk[3]
        xk = hitk[0]
        yk =hitk[1]

        radiationLength=kp.LrXe/kp.Pr

        zslice=KFZSlice(radiationLength,zki,zkf,edepk)
        hit =KFVector([xk,yk])
        cov=KFMatrix([[kp.sigma_x**2,0],
                   [0,kp.sigma_y**2]])
                
        measurement = KFMeasurement(hit,cov)
        node =KFNode(measurement,zslice)
        Nodes.append(node)
        Dz.append(zkf-zki)
        Zi.append(zki)

    if debug >= Debug.draw.value:
        mpldv = MPLDrawVector(Dz)
        mpldv.Draw()
        mpldv = MPLDrawVector(Zi)
        mpldv.Draw()
    setup = KFSetup(Nodes) 
    system =KFSystem(setup)
    return system

def KFLinealInitState(hit0,betaMax,smearVector):
        """
        Initial State for fitter KFLineal 
        """

        x0 = hit0[0]
        y0 = hit0[1]
        ux0,uy0 = betaMax.DirectionTangents()
        
        sx0 = smearVector[0]
        sy0 = smearVector[1]
        sux0 =0.5*ux0 # arbitraryly large relative value
        suy0 =0.5*uy0

        x0x0 = sx0**2
        y0y0 = sy0**2     
        ux0ux0 =sux0**2
        uy0uy0 =suy0**2
        x0ux0 = sx0*sux0
        x0uy0 = sx0*suy0
        y0ux0 = sy0*sux0
        y0uy0 = sy0*suy0
        x0y0 = sx0*sy0
        ux0uy0 = sux0*suy0

        traj0 = KFVector([x0,y0,ux0,uy0])
        cov0 = KFMatrix([
                    [x0x0,x0y0,x0ux0,x0uy0],
                    [x0y0,y0y0,y0ux0,y0uy0],
                    [x0ux0,y0ux0,ux0ux0,ux0uy0],
                    [x0uy0,y0uy0,ux0uy0,uy0uy0]
                    ])


        return (traj0,cov0)    
def DrawHits(betaMax):
    """
    Given a KMCParticle, draw hits, using root or Matplotlib
    """

    print "drawing"
    trueHits = betaMax.TrueHits()
    smearHits = betaMax.SmearedHits()
    scale =.6*kp.Pr
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
        drawHits.DrawMeasurements(draw='2D')
        drawHits.DrawMeasurements(draw='3D')

        drawHits.DrawHits(draw='2D')
        drawHits.DrawHits(draw='3D')

        drawHits.DrawAll(draw='2D')
        

