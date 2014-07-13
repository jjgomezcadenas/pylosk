
from KFBase import KFVector, KFMatrix
from KFMeasurement import KFMeasurement 
from KFZSlice import KFZSlice
from KFNode import KFNode 
from KFState import KFState
from KFSetup import KFSetup 
from KFSystem import KFSystem
from KMCParticle import KMCParticle

import sys

from math import *
from abc import ABCMeta, abstractmethod

from LogKalmanFilter import *

"""
Kalman Filter is the main module that perfoms the Kalman Filter fit. 
JJGC, May, 2014
"""


class KalmanFilter(object):
    __metaclass__ = ABCMeta
    """
    Represents a KF fitter
    
    Takes a System and additional conditions as extra arguments   
    P0 is the initial momentum (needed for the computation of MS) 
    """
    
    def __init__(self,name,P0=2.9):
        """
        Inits the fitter with a name and the initial
        momentum of the track (or a guess)
        """

        self.name = name
        self.P0 = P0
    

    def SetSystem(self,system):
        """
        Sets the system
        """

        if isinstance(system, KFSystem) != True:
            print "error, must set with an instance of KFSystem " 
            sys.exit(-1)

        self.system = system

    def GetSystem(self):
        """
        Returns the System
        """
        return self.system

    def GetStates(self):
        """
        Returns the States of the system
        """
        return self.system.States

    def GetNodes(self):
        """
        Returns the States of the system
        """
        return self.system.Nodes

    def SetKMCParticle(self,kmcparticle):
        """
        Sets the KMCParticle
        """

        if isinstance(kmcparticle, KMCParticle) != True:
            print "error, must set with an instance of KMCParticle " 
            sys.exit(-1)

        self.kmcparticle = kmcparticle

    def GetKMCParticle(self):
        """
        Returns the KMCParticle
        """
        return self.kmcparticle

    def GetInitialMomentum(self):
        """
        Returns P0
        """
        return self.P0
    
    def FitterName(self):
        """
        Provides a string identifying the type of fitter
        """
        return self.name

    @abstractmethod
    def H(self,k):
        """
        H matrix, to be provided by specific implementation of fitter
        (straight line, helix)
        """
        return

    @abstractmethod
    def TransportMatrix(self,k):
        """
        Transport matrix Fk to be provided by specific 
        implementation of fitter
        """
        return 
    @abstractmethod   
    def MultipleScatteringMatrix(self,k):
        """
        The multiple scattering matrix Qk to be provided 
        by specific implementation of fitter
        """
        return 

    @abstractmethod
    def SetInitialState(self):
        """
        Set the initial state
        """   
        return
        
    @abstractmethod
    def SetupFilter(self,k):
        """
        Return a state at point k that is "set up" using the measurements
        and the current momentum of the KalmanFilter.  Note the current
        momentum is the value of self.P, which changes as the filter is
        applied.
        """
        return


    def Predict(self,k=1):
        """
        Prediction from k-1:

            ak,k-1 = FK-1*ak-1
            Ck,k-1 = Fk-1 Ck-1 Fk-1^T + Qk-1

            where (ak-1,Ck-1) represent the state in node k-1
            notice:

            formula does not apply for k=0 (no prediction for initial site)
            ak,k-1 = aP (P = predicted)
            ak-1 = aF (F = Filtered)
            Ck,k-1 = CP (C predicted)
            Ck-1 = CF (C filtered from previous state)
            Fk-1 : transport from k-1 to k (evaluated from dz) =F
            Qk-1 : MS matrix, evaluated from dz

            aP = F*aF
            CP = F*CF*F^T + Q

        """

        lgx.debug("-->Predict: k = {0}".format(k))
            

        if k== 0:
            print "Can't predict at site 0 k must be >=1"
            sys.exit(-1)

        stateDict=self.system.States[k-1]
        state =stateDict["F"]

        aF = state.V
        CF = state.Cov
        F = self.TransportMatrix(k)   # the method should be constructed to give F(k-1)
        Q =self.MultipleScatteringMatrix(k) # the method should be constructed to give Q(k-1)
        FT = F.Transpose()
        aP = F*aF
        CPt = F*CF*FT 
        CP = F*CF*FT + Q

        lgx.debug("Filtered from previous state->\n aF={0}\n CF={1}".
            format(aF,CF))

        lgx.debug("Transport and MS->\n F={0} \n FT={1} \n Q={2}".
            format(F,FT,Q))

        lgx.debug("Predict-> \n aP=F*aF->{0}\n MS matrix = {1}".
            format(aP,Q))

        lgx.debug("Predict-> CPt= F*CF*FT->{0} CP= F*CF*FT+Q ->{1}".
            format(CPt,CP))

        pstate = KFState(aP,CP)
        pstate.Chi2=1e+6
        self.system.States.append({"P":pstate,"F":0,"S":0})

        cond_pause(debug)

        return 0

    def Filter(self,k=1):
        """
        Filter at k:

        Ck = [(Ck,k-1)^-1 + H^T Gk H]^-1

        where Ck is the filtered covariance matrix obtained from:
        Ck,k-1: the predicted covariance matrix
        H: the matrix that transform trajectory to measurement: mk = H*ak
        G: the inverse of the covariance of measurements: G = V^-1

        then one computes the Kalman Gain Matrix:

        K = Ck*H^T*G

        and the residual:

        R = mk - H * ak,k-1

        where mk is the vector of measurements (hit), ak,k+1 is the predicted trajectory
        and H * ak,k-1 gives the hit predicted.

        Finally, the filtered state is:

        ak = ak,k-1 + K * R

        We call:
        ak = aF  
        ak,k-1 = aP
        Ck = CF
        Ck^-1 = CFI 
        Ck,k-1 = CP
        Ck, k-1^-1 =CPI 
        """

        lgx.debug("-->Filter : k = {0}".format(k))

        if k== 0:
            print "Can't filter at site 0 k must be >=1"
            sys.exit(-1)
        
        stateDict=self.system.States[k]
        state =stateDict["P"]  #predicted
        aP = state.V  # aP = ak,k+1
        CP = state.Cov

        lgx.debug("Predicted State->  aP ={0} \n  CP ={1}".
            format(aP,CP))

        m = self.system.Nodes[k].Measurement.V
        G = self.system.Nodes[k].Measurement.Cov.Inverse()

        lgx.debug("Measurement-> \n  m ={0} \n G = V^-1 ={1}".
            format(m,G))

        cond_pause(debug)

        CPI = CP.Inverse()
        H = self.H(k)
        HT =H.Transpose()

        lgx.debug("Matrices-> \n Ck,k-1^-1 =CPI ={0} \n H ={1} \n HT ={2}".
                 format(CPI,H,HT))

        cond_pause(debug)

        CFI = CPI + HT*G*H  # CFI = Ck^-1
        CF = CFI.Inverse() #CF = Ck : filtered
        K = CF*HT*G  # G= V^-1: K : Gain Matrix

        r2 = H*aP   # transforms predicted state

        lgx.debug("Filt Cov-> \n Ck,^-1 =CFI ->{0} \n Ck = CF = {1} \n K =CF*HT*G ->{2}".
                 format(CFI,CF,K))

        cond_pause(debug)

        lgx.debug("-->r2 =H*aP ={0}".
                 format(r2))
        
        r = m - H*aP # residual wrt predicted

        lgx.debug("-->r = m - H*aP ={0}".
                 format(r))

        cond_pause(debug)

        r2 = K*r # K =(4x2) r = (2x1): KfMatrix*KFvector -->(KFvector) (4x1)

        lgx.debug(" r2 =K*r ={0}  ".
                 format(r2))

        cond_pause(debug)
        
        aF = aP + K* r # aF = ak : filtered

        lgx.debug("aF = aP + K* r-> \n aP ={0} \n aF ={1} \n ".
                 format(aP,aF))

        cond_pause(debug)

        r2 = H*aF

        lgx.debug("r2 =H*aF ={0} \n".
                 format(r2))

        R = m - H*aF # residual wrt Filtered

        lgx.debug("R = m - H*aF-> \n R ={0}  \n ".
                 format(R))

        RT = R.Transpose()
        
        lgx.debug("RT  ={0}, G={1}  \n ".
                 format(RT,G))

        cond_pause(debug)
        
        Rf = aF - aP          # residual (filtered - predicted)
        RfT = Rf.Transpose()
        
        chi2meas = RT*G*R     # chi2 contribution due to measurement
        chi2ms = RfT*CPI*Rf   # chi2 contribution due to MS
        
        chi2 = (chi2meas[0] + chi2ms[0])/2.  # chi2/2. is chi2/dof

        lgx.debug("Residuals-> \n r ={0} \n R = {1} ".
                 format(r,R))

        lgx.info("Chi2Meas ={0} Chi2Ms = {1}".
                 format(chi2meas,chi2ms))

        lgx.info("Filtered State-> \n aF ={0} \n Chi2 = {1} ".
                 format(aF,chi2))

        # Construct the filtered state.
        fstate = KFState(aF,CF)
        fstate.Chi2=chi2
        stateDict["F"] = fstate

        # Save the MS part of the chi2 in the predicted state.
        stateDict["P"].Chi2 = chi2ms[0]/2; 

        cond_pause(debug)
        
        return 0


    def SigmaThetaMs(self,P,L):
        """
        sigma(theta_ms) = 13.6 (Mev)/(beta*P)*Sqrt(L)*(1+0-038*log(L))
        L in radiation length 
        """
        beta = self.Beta(P)
        if(beta <= 0): return 0.;
 
        try:
            tms = ((13.6)/(P*1.*beta))*sqrt(L*1.)*(1 + 0.038*log(L*1.))
        except ValueError:
            lgx.debug("+++KalmanFilter:SigmaThetaMs+++: Error calculating tms: P = {0}, beta = {1}, L = {2}".format(P,beta,L))
            raise;

        lgx.info("SigmaThetaMs->  L={0} P={1} beta={2}, tms={3}".
            format(L,P,beta,tms))

        cond_pause(debug)

        return tms


    def Sigma2ThetaMs(self,P,L):
        """
        sigma2(theta_ms) 
        """
        tms =self.SigmaThetaMs(P,L)

        return tms**2

    def Beta(self, P):
        """
        beta = P/E
        """
        E = sqrt(P**2+0.511**2)
        beta = P/E

        # lgx.debug("Beta-> P ={0}, E={1} beta={2}".
        #     format(P,E,beta))
        return beta

    def CorrectP(self, P0, edep):
        E0 = sqrt(P0**2+0.511**2)
        E = E0-edep
        if(E > 0.511):
            P= sqrt(E**2-0.511**2)
        else:
            P= 0.

        return P



# class KFLineal(KalmanFilter):
#     """
#     Represents a KF fitter for the case of a straight line    
#     """
        
#     def __init__(self,name,P0=2.9):
#         """
#         """

#         KalmanFilter.__init__(self,name)
        
#         self.P = P0 # in MeV
#         lgx.info("KFLineal__init__ -> P0 ={0}".
#                format(P0))


#     def H(self,k):
#         """
#         H matrix
#         """
#         return KFMatrix([[1.,0.,0.,0.],[0.,1.,0.,0.]])


#     def SetInitialState(self):
#         """
#         Set the initial state
#         """   
        
#         p0 = self.P;
#         v0 = self.system.Nodes[k].Measurement.V;
#         v1 = self.system.Nodes[k+1].Measurement.V;
#         zi0 = self.system.Nodes[k].ZSlice.Zi
#         zi1 = self.system.Nodes[k+1].ZSlice.Zi
#         zf0 = self.system.Nodes[k].ZSlice.Zf
#         zf1 = self.system.Nodes[k+1].ZSlice.Zf
        
#         # Set the initial trajectory (get initial tangent angles from first two points).
#         x0 = v0[0]; x1 = v1[1]
#         y0 = v0[0]; y1 = v1[1];
#         z0 = (zi0 + zf0)/2.; z1 = (zi1 + zf1)/2.;
#         z02 = z0**2;
#         mag = sqrt((x1 - x0)**2 + (y1-y0)**2 + (z1-z0)**2);
#         ux = (x1 - x0)/mag;
#         uy = (y1 - y0)/mag;
#         uz = (z1 - z0)/mag;
#         tanX0 = ux/uz; tanY0 = uy/uz;
    
#         p1 = trk_x0[0] - z0*tanX0;
#         p2 = trk_y0[0] - z0*tanY0;
#         p3 = tanX0;
#         p4 = tanY0;
#         traj0 = KFVector([p1,p2,p3,p4]);

#         x0 = hit0[0]
#         y0 = hit0[1]
#         ux0,uy0 = betaMax.DirectionTangents()
        
#         sx0 = smearVector[0]
#         sy0 = smearVector[1]
#         sux0 =0.5*ux0 # arbitraryly large relative value
#         suy0 =0.5*uy0

#         x0x0 = sx0**2
#         y0y0 = sy0**2     
#         ux0ux0 =sux0**2
#         uy0uy0 =suy0**2
#         x0ux0 = sx0*sux0
#         x0uy0 = sx0*suy0
#         y0ux0 = sy0*sux0
#         y0uy0 = sy0*suy0
#         x0y0 = sx0*sy0
#         ux0uy0 = sux0*suy0

#         traj0 = KFVector([x0,y0,ux0,uy0])
#         cov0 = KFMatrix([
#                     [x0x0,x0y0,x0ux0,x0uy0],
#                     [x0y0,y0y0,y0ux0,y0uy0],
#                     [x0ux0,y0ux0,ux0ux0,ux0uy0],
#                     [x0uy0,y0uy0,ux0uy0,uy0uy0]
#                     ])


         
#         state0 = KFState(traj0,cov0)

#         self.system.States.append({"P":state0,"F":state0,"S":0})
#     def TransportMatrix(self,k):
#         """
#         The Transport Matrix F allows to transport the State from site i to j
#         """

#         dz =self.system.Nodes[k].ZSlice.dz 
        
#         return KFMatrix([   [1.,0,dz,0.],
#                             [0,1.,0.,dz],
#                             [0,0,1.,0],
#                             [0,0,0,1.]]) 

#     def MultipleScatteringMatrix(self,k):
#         """
#         The multiple scattering matrix Qk
#         """
#         dz =self.system.Nodes[k].ZSlice.dz
#         dz2=dz*dz
#         edep = self.system.Nodes[k].ZSlice.Edep*1.
#         Lr = self.system.Nodes[k].ZSlice.Lr
#         L=abs(dz)/Lr

#         stateDict=self.system.States[k-1]
#         state =stateDict["F"]

#         ak = state.V
#         p1 = ak[0]
#         p2 = ak[1]
#         p3 = ak[2]
#         p4 = ak[3]

#         lgx.debug("SlineKF:MultipleScatteringMatrix -> ak ={0}".
#             format(ak))

#         lgx.debug(" -> dz ={0} dz2 = {1} edep ={2}".
#             format(dz,dz2,edep))

#         lgx.debug(" -> Lr ={0} L = {1} P ={2}".
#             format(Lr,L,self.P))

#         self.P =self.CorrectP(self.P, edep) 

#         lgx.debug(" -> P (after correct) ={0}".
#             format(self.P))
        
#         s2tms =self.Sigma2ThetaMs(self.P,L)
#         lgx.debug(" -> s2tms ={0}".
#             format(s2tms))

#         cond_pause(debug)

#         p3p3 =s2tms*(1+p3**2)*(1+p3**2+p4**2)
#         p4p4 = s2tms*(1+p4**2)*(1+p3**2+p4**2)
#         p3p4=s2tms*p3*p4*(1+p3**2+p4**2)

#         lgx.debug(" -> p3p3 ={0} p4p4 ={1} p3p4 ={2} ".
#             format(p3p3,p4p4,p3p4))

#         cond_pause(debug)

#         covList = [
#         [dz2*p3p3,dz2*p3p4,-dz*p3p3,-dz*p3p4],
#         [dz2*p3p4,dz2*p4p4,-dz*p3p4,-dz*p4p4],
#         [-dz*p3p3,-dz*p3p4,p3p3,p3p4],
#         [-dz*p3p4,-dz*p4p4,p3p4,p4p4]
#                    ]

#         # covList = [
#         # [0,0,0,0],
#         # [0,0,0,0],
#         # [0,0,p3p3,p3p4],
#         # [0,0,p3p4,p4p4]
#         #            ]
        

#         cov=KFMatrix(covList)

#         s = """

#         cov=KFMatrix([
#         [dz2*p3p3 = {0},dz2*p3p4 = {1},-dz*p3p3 = {2},-dz*p3p4 = {3}],
#         [dz2*p3p4 = {4},dz2*p4p4 = {5},-dz*p3p4 = {6},-dz*p4p4 = {7}],
#         [-dz*p3p3 = {8},-dz*p3p4 = {9},p3p3 = {10},p3p4 = {10}],
#         [-dz*p3p4 = {12},-dz*p4p4 = {13},p3p4 = {14},p4p4 = {15}]
#                    ])
#         """
#         lgx.debug(s.format(dz2*p3p3,dz2*p3p4,-dz*p3p3,-dz*p3p4,
#         dz2*p3p4,dz2*p4p4,-dz*p3p4,-dz*p4p4,
#         -dz*p3p3,-dz*p3p4,p3p3,p3p4,
#         -dz*p3p4,-dz*p4p4,p3p4,p4p4))

#         return cov
    

# def Setup():
#     LrXe = 15300. # radiation length of xenon in mm
#     Pr = 10.      #pressure in bar

#     Nodes=[]
#     for i in range(1,6):
#         zslice=KFZSlice(LrXe/Pr,i,i+1,0.1)
    
#         x0=5+i
#         y0=15+i
#         hit =KFVector([x0,y0])
#         sigma_x=1.0
#         sigma_y=1.0
#         cov=KFMatrix([[sigma_x**2,0],
#                    [0,sigma_y**2]])

#         measurement = KFMeasurement(hit,cov)

#         node =KFNode(measurement,zslice)
#         Nodes.append(node)

#     setup = KFSetup(Nodes) 
#     return setup

# def InitialState():
#     x0 = 5.
#     y0 = 10.
#     ux0 = 0.25
#     uy0 = 0.25
#     sx0 = 0.25
#     sy0 = 0.25
#     x0x0 = sx0**2
#     y0y0 = sy0**2

#     sux0 =1.
#     suy0 =1.
#     ux0ux0 =sux0**2
#     uy0uy0 =suy0**2

#     x0ux0 = sx0*sux0
#     x0uy0 = sx0*suy0

#     y0ux0 = sy0*sux0
#     y0uy0 = sy0*suy0

#     x0y0 = sx0*sy0
#     ux0uy0 = sux0*suy0

#     traj0 = KFVector([x0,y0,ux0,uy0])


#     cov0 = KFMatrix([
#                     [x0x0,x0y0,x0ux0,x0uy0],
#                     [x0y0,y0y0,y0ux0,y0uy0],
#                     [x0ux0,y0ux0,ux0ux0,ux0uy0],
#                     [x0uy0,y0uy0,ux0uy0,uy0uy0]
#                     ])

#     return (traj0,cov0) 
    
if __name__ == '__main__':
    print "OK"

#     setup = Setup()  #creates a dummy setup
#     print "Setup =",setup

#     system =KFSystem(setup)  #System takes a setup

#     slkf = KFLineal("KFLineal",P0=2.9) 
#     slkf.SetSystem(system)

#     a0,C0 = InitialState() #initial state

#     state0 = KFState(a0,C0)
#     print "Initial State =",state0

#     slkf.SetInitialState(state0) #Init fitter with initial state
    
#     print "slkf States after init"
#     print slkf.System.States

#     slkf.Predict(k=1) #predicts from state 0 (init) to state 1

#     print "slkf States after predict"
#     print slkf.System.States

#     slkf.Filter(k=1) #Filters prediction in state 1

#     print "slkf States after filter"
#     print slkf.System.States
