"""
KFWolinFilter.py

A KalmanFilter following the state conventions described in Wodin et. al.

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

from LogKFWolinFilter import *

class KFWolinFilter(KalmanFilter):
    """
    A KF filter following Wodin et. al.  
    """
        
    def __init__(self,name,P0=2.393):
        """
        inits the fitter
        """

        KalmanFilter.__init__(self,name)
        
        self.P = P0 # in MeV
        lgx.info("KFWolinFilter__init__ -> P0 ={0}".
            format(P0))
    
    def SetupFilter(self,k):
        """
        Constructs a state with momentum equal to the current
        momentum of the Kalman filter and using the directional fit
        from slices k to k+4.
        """
        
        p0 = self.P
        nodes = self.GetNodes()
        nfpts = 4
        
        # Set up the lists.
        fpt_x = []; fpt_y = []; fpt_z = []
        for tk in range(k,k+nfpts):
            
            if(tk < len(nodes)):
                fV = nodes[tk].Measurement.V
                fx = fV[0]; fy = fV[1]
                fzi = nodes[tk].ZSlice.Zi
                fzf = nodes[tk].ZSlice.Zf
                fz = (fzi + fzf)/2.
                
                fpt_x.append(fx)
                fpt_y.append(fy)
                fpt_z.append(fz)
                
        # Set the z0 and z0**2 for (zi of the slice) for calculation of MS.
        z0 = nodes[k].ZSlice.Zi
        z02 = z0**2

        # Fit the points to a line to predict the next direction vector.
        # Direction vector of line is vv[0] returned by np.linalg.svd:
        #   code from: http://stackoverflow.com/questions/2298390/fitting-a-line-in-3d
        xvals = np.array(fpt_x)
        yvals = np.array(fpt_y)
        zvals = np.array(fpt_z)
        dmatrix = np.concatenate((xvals[:, np.newaxis], 
                                  yvals[:, np.newaxis], 
                                  zvals[:, np.newaxis]), axis=1)
        dmean = np.array([fpt_x[0],fpt_y[0],fpt_z[0]]) # ensure the line passes through the point k
        uu, dd, vv = np.linalg.svd(dmatrix - dmean)
        uvec = vv[0]

        # Calculate the length of the next step.
        Lr = nodes[k].ZSlice.Lr
        dz = nodes[k].ZSlice.dz
        L = abs(dz)/Lr
        
        # Set the initial trajectory (get initial tangent angles from first two points).
        x0 = fpt_x[0]
        y0 = fpt_y[0]
        ux = uvec[0]
        uy = uvec[1]
        uz = uvec[2]
        
        tanX0 = ux/uz; tanY0 = uy/uz
        p1 = x0 - z0*tanX0
        p2 = y0 - z0*tanY0
        p3 = tanX0
        p4 = tanY0
        traj0 = KFVector([p1,p2,p3,p4])

        # Compute the initial MS covariance matrix.
        s2tms = self.SigmaThetaMs(p0,L)
        s2tms = s2tms**2
        p3p3 = s2tms*(1+p3**2)*(1+p3**2+p4**2)
        p4p4 = s2tms*(1+p4**2)*(1+p3**2+p4**2)
        p3p4= s2tms*p3*p4*(1+p3**2+p4**2)

        cov0 = KFMatrix([
            [z02*p3p3,z02*p3p4,-z0*p3p3,-z0*p3p4],
            [z02*p3p4,z02*p4p4,-z0*p3p4,-z0*p4p4],
            [-z0*p3p3,-z0*p3p4,p3p3,p3p4],
            [-z0*p3p4,-z0*p4p4,p3p4,p4p4]
            ])

        # Create the state object.
        state0 = KFState(traj0,cov0)
        
        return state0

    def SetInitialState(self):
        """
        Set the initial state
        """   
        state0 = self.SetupFilter(0)
        self.GetStates().append({"P":state0,"F":state0,"S":0})
    
    def H(self,k):
        """
        H matrix at state k
        """
        dz =self.system.Nodes[k].ZSlice.dz;
        z0 =self.system.Nodes[k].ZSlice.Zi + dz/2.;
        return KFMatrix([[1.,0.,z0,0.],[0.,1.,0.,z0]])


    def TransportMatrix(self,k):
        """
        The Transport Matrix F allows to transport the State from site i to j
        """
        
        return KFMatrix([[1.,0.,0.,0. ],
                            [0.,1.,0.,0.],
                            [0.,0.,1.,0.],
                            [0.,0.,0.,1]]) 

    def MultipleScatteringMatrix(self,k):
        """
        The multiple scattering matrix Qk
        """
        dz =self.system.Nodes[k].ZSlice.dz
        z0 =self.system.Nodes[k].ZSlice.Zi
        z02=z0*z0;
        edep = self.system.Nodes[k].ZSlice.Edep*1.
        Lr = self.system.Nodes[k].ZSlice.Lr
        L=abs(dz)/Lr

        stateDict=self.system.States[k-1]
        state =stateDict["F"]

        ak = state.V
        p1 = ak[0]
        p2 = ak[1]
        p3 = ak[2]
        p4 = ak[3]

        lgx.debug("WodinKFilter:MultipleScatteringMatrix -> ak ={0} type ={1}".
            format(ak,type(ak)))

        lgx.debug(" -> p1 ={0} p2 ={1} p3 ={2} p4 ={3}, type(p1)={4}".
            format(p1,p2,p3,p4,type(p1)))

        lgx.debug("WodinKFilter:MultipleScatteringMatrix -> dz ={0} dz2 = {1} edep ={2}".
            format(dz,dz*dz,edep))

        lgx.debug("WodinKFilter:MultipleScatteringMatrix -> Lr ={0} L = {1} P ={2}".
            format(Lr,L,self.P))

        self.P =self.CorrectP(self.P, edep)
        
        lgx.debug("WodinKFilter:MultipleScatteringMatrix -> P (after correct) ={0}".
            format(self.P))
        
        s2tms =self.Sigma2ThetaMs(self.P,L)

        lgx.debug("WodinKFilter:MultipleScatteringMatrix -> s2tms ={0} type ={1}".
            format(s2tms,type(s2tms)))

        p3p3 =s2tms*(1+p3**2)*(1+p3**2+p4**2)
        p4p4 = s2tms*(1+p4**2)*(1+p3**2+p4**2)
        p3p4=s2tms*p3*p4*(1+p3**2+p4**2)

        cov=KFMatrix([
        [z02*p3p3,z02*p3p4,-z0*p3p3,-z0*p3p4],
        [z02*p3p4,z02*p4p4,-z0*p3p4,-z0*p4p4],
        [-z0*p3p3,-z0*p3p4,p3p3,p3p4],
        [-z0*p3p4,-z0*p4p4,p3p4,p4p4]
                   ])

        return cov
    

def Setup():
    LrXe = 15300. # radiation length of xenon in mm
    Pr = 10.      #pressure in bar

    Nodes=[]
    for i in range(1,6):
        zslice=KFZSlice(LrXe/Pr,i,i+1,0.1)
    
        x0=5+i
        y0=15+i
        hit =KFVector([x0,y0])
        sigma_x=0.1
        sigma_y=0.2
        cov=KFMatrix([[sigma_x**2,0],
                   [0,sigma_y**2]])

        measurement = KFMeasurement(hit,cov)

        node =KFNode(measurement,zslice)
        Nodes.append(node)

    setup = KFSetup(Nodes) 
    return setup

def InitialState():
    x0 = 5.
    y0 = 10.
    ux0 = 0.25
    uy0 = 0.25
    sx0 = 0.1
    sy0 = 0.1
    x0x0 = sx0**2
    y0y0 = sy0**2

    sux0 =0.1
    suy0 =0.1
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
    
if __name__ == '__main__':

    setup = Setup()  #creates a dummy setup
    print "Setup =",setup

    system =KFSystem(setup)  #System takes a setup

    slkf = KFWolinFilter("KFWolin",P0=2.9) #KF Fitter 
    slkf.SetSystem(system)
    slkf.SetInitialState()
    # a0,C0 = InitialState() #initial state

    # state0 = KFState(a0,C0)
    # print "Initial State =",state0

    #slkf.SetInitialState(state0) #Init fitter with initial state
    
    print "slkf States after init"
    print slkf.GetStates()

    slkf.Predict(k=1) #predicts from state 0 (init) to state 1

    print "slkf States after predict"
    print slkf.GetStates()

    slkf.Filter(k=1) #Filters prediction in state 1

    print "slkf States after filter"
    print slkf.GetStates()
