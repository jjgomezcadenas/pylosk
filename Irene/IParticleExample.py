from ROOT import *
import IParam as ip
import numpy as np
import sys,getopt,os

from IParticle import IParticle
from TDrawHits import TDrawHits
from MPLDrawHits import MPLDrawHits

import logging 
logging.basicConfig(level=logging.DEBUG)



def cond_pause(debug_level):
    if debug_level <= ip.debug:
        s=raw_input("return to continue")

def Analyze(ievt,n_event):
    ipmax,itrkmax = SelectEMax(ievt)

    smearVector =np.array([ip.sigma_x, ip.sigma_y, 0.]) 
    
    scale = ip.Pr/2.
    betaMax = IParticle(ipmax,itrkmax,smearVector,ip.sample,Pressure=ip.Pr)

    print """
    betaMax ={0}
    """.format(betaMax)
    cond_pause(ip.Debug.chat.value)


    trueHits = betaMax.TrueHits()
    smearHits = betaMax.SmearedHits()
    lowerBin =np.array([-700.,-700.,-700.])/scale
    upperBin =np.array([700.,700.,700.])/scale

    if ip.DrawRoot == True:
        drawHits = TDrawHits((lowerBin,upperBin),trueHits,smearHits)
        drawHits.DrawMeasurements(draw='2D')
        drawHits.DrawMeasurements(draw='3D')

        drawHits.DrawHits(draw='2D')
        drawHits.DrawHits(draw='3D')

        drawHits.DrawAll(draw='2D')
        drawHits = MPLDrawHits((lowerBin,upperBin),trueHits,smearHits)

    if ip.DrawMPL == True:
        drawHits.DrawMeasurements(draw='2D')
        drawHits.DrawMeasurements(draw='3D')

        drawHits.DrawHits(draw='2D')
        drawHits.DrawHits(draw='3D')

        drawHits.DrawAll(draw='2D')
    

def main(argv):
    
    inputFile=inputFile
    inputDir=ip.inputDir
    n_events=ip.n_events
    
    print """
        Input file ={0}
        Input dir ={1}
        Events to run ={2}
        Path to data ={3}
        Debug Level ={4}
        """.format(inputFile,inputDir,n_events,os.environ['DATA'],ip.debug)
        
        
    source=os.environ['DATA']+'/'+inputFile

    print "Opening file =%s and loading TTree"%source
        
    fFile = TFile.Open(source)
    fEvtTree = fFile.Get(ip.gsEvtTree)
    n_entries = fEvtTree.GetEntries()

   
    print "number of entries in file =%d"%(n_entries)

    timer = TStopwatch()
    timer.Start()	

    print "Loading Irene Event"

    ievt =irene.Event()  #create an irene event
    fEvtTree.SetBranchAddress(ip.gsEvtBranch, ievt)
    nb = 0 #number of bytes read

    n_run = min(n_entries,int(n_events))

    print "number of events to run = %d"%n_run
    s=raw_input("return to continue")

    
    for n_event in range(0,n_run):
        print "Reading event = %d"%n_event
        s=raw_input("return to continue")
            
        nb += fEvtTree.GetEntry(n_event)
        Analyze(ievt,n_event)
        

    timer.Stop() 
    print "I am out of the event loop"
    mbytes = 0.000001*nb #in megabytes
    rtime = timer.RealTime()
    ctime = timer.CpuTime()
    print "RealTime=%5.3f seconds, CpuTime=%5.3f seconds\n"%(rtime,ctime)
    print "Read %5.2f Mbytes/Realtime seconds\n"%(mbytes/rtime)
    print "Read %5.2f Mbytes/Cputime seconds\n"%(mbytes/ctime)
    print "%d events and %d bytes read.\n"%(n_event,nb);
    fFile.Close();

    
    s=raw_input("return to continue") 
    print "done"

def SelectEMax(ievt):
    """
        Selects the electron with max energy (betaMax)
    """
    
    itrks= ievt.GetTracks()
    n_itrks = itrks.GetEntries()
  
    ipmax = 0
    itrkmax=0
    for it in range(0,n_itrks):
        itrk = itrks.At(it) 
        ipart= itrk.GetParticle()

        if ipart.IsPrimary() == False:
            continue

        T= (ipart.Energy()-ipart.GetMass())*1e+3
        if abs(T-ip.TMAX)< ip.eps:
            ipmax = ipart
            itrkmax=itrk

    return ipmax,itrkmax

def Setup(Hits):
    """
    Creates a KF Setup 
    """

    logging.info("Now in Method: Setup")
    Nodes=[]
    
    for k in range(0,len(Hits)):
        hitk =Hits[k]
        hitkp1 =Hits[k+1]
        zki =hitk[2]
        zkf = hitkp1[2]
        edepk = hitk[3]
        xk = hitk[0]
        yk =hitk[1]
        zslice=KFZSlice(ip.LrXe/ip.Pr,zki,zkf,edepk)
        hit =KFVector([xk,yk])
        cov=KFMatrix([[ip.sigma_x**2,0],
                   [0,ip.sigma_y**2]])
                
        measurement = KFMeasurement(hit,cov)

        node =KFNode(measurement,zslice)

        
        Nodes.append(node)

    logging.info("Setup: number of nodes ={0}".format(len(Nodes)))
    
    setup = KFSetup(Nodes) 
    return setup

if __name__ == '__main__':
    gSystem.Load("libirene")
    main(sys.argv[1:])


