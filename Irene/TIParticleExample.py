from ROOT import *

import IParam as ip
import numpy as np
import sys,getopt,os

from TIParticle import *

def Analyze(ievt,n_event):
    ipmax,itrkmax = SelectEMax(ievt)

    smearVector =np.array([1.0, 1.0, 0.]) 
    sample=5
    Pr =5
    betaMax = TIParticle(ipmax,itrkmax,smearVector,sample,Pressure=Pr)

    print """
    betaMax ={0}
    """.format(betaMax)
    s=raw_input("return to continue") 

    betaMax.DrawHits(draw='2D',view='True')
    betaMax.DrawHits(draw='3D',view='True')
    betaMax.DrawHits(draw='2D',view='Smear')
    betaMax.DrawHits(draw='3D',view='Smear')
    betaMax.DrawHits(draw='2D',view='Both')
    betaMax.DrawHits(draw='3D',view='Both')

    #Get true hits
    trueHits = betaMax.TrueHits()
    
    

def main(argv):
    
    inputFile="MagBox_Xe_5atm_00tesla.e2447.next"
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


if __name__ == '__main__':
    gSystem.Load("libirene")
    main(sys.argv[1:])


