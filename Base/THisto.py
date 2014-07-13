"""
Module THisto 
"""
import numpy as np
import sys
from ROOT import TH1F, TH2F,TH3F
from ROOT import TCanvas, TPad, TPaveLabel, TPaveText, TColor
from ROOT import gROOT, gStyle, gRandom
gStyle.SetOptStat(0);
gStyle.SetPalette(1);
gStyle.SetCanvasColor(33);
gStyle.SetFrameFillColor(18);


class THisto(object):
    """
    An root-based histogram server 
    """

    def __init__(self):
            
        """
        H1 is a dictionary of 1-dim histos
        """
        self.H1={}
        self.H2={}

    def BookH1(self,name,title,nbin,xmin,xmax):
        self.H1[name]=TH1F(name,title,nbin,xmin,xmax)

    def BookH2(self,name,title,nbinx,xmin,xmax,nbiny,ymin,ymax):
        self.H2[name]=TH2F(name,title,nbinx,xmin,xmax,nbiny,ymin,ymax)
        
    def FillH1(self,name,x,w=1):
        self.H1[name].Fill(x,w)

    def FillH2(self,name,x,y,w=1):
        self.H2[name].Fill(x,y,w)

    def __str__(self):
        s="H1 ="+self.H1.__str__()
        s+="\nH2="+self.H2.__str__()
        return s

    def __repr__(self):
        return self.__str__()


    def DrawAll(self):
        """
        Draws histos 
        """
        
        Histos=[]

        for name in self.H1.keys():
            Histos.append(self.H1[name])

        for name in self.H2.keys():
            Histos.append(self.H2[name])

        self.__DrawHistos(Histos,1,1)


    def DrawList(self,histoList,xd=1,yd=1):
        """
        Draws histos 
        """

        Histos=[]
        for name in histoList:
            if name in self.H1:
                Histos.append(self.H1[name])
            elif name in self.H2:
                Histos.append(self.H2[name])
            else:
                print "histo name ={0} not found".format(name)
                sys.exit(-1)

        self.__DrawHistos(Histos,xd,yd)


    def __DrawHistos(self,Histos,xd,yd):

        c1 = TCanvas( 'c1', 'Histograms', 200, 10, 600, 800 )

        l = len(Histos)
        div = xd*yd

        if not float(div)%float(l)==0:
            yd = int(l*1./xd*1.)

        
        c1.Divide(xd,yd)
        
        for i in range (1,l+1):
            c1.cd(i)
            Histos[i-1].Draw()
        

        s=raw_input("return to continue")


if __name__ == '__main__':
    gRandom.SetSeed()

    th = THisto()
    th.BookH1("s1","This is the first signal",100,-4,4)
    th.BookH1("s2","This is the second signal",100,-4,4)
    th.BookH2("s1s2","This is the first versus the econd signal",100,-4,4,100,-4,4)

    for i in range(0,500):
        xs1   = gRandom.Gaus(-0.5,0.5)
        xs2   = gRandom.Landau(1,0.15)
        th.FillH1("s1",xs1,0.3)
        th.FillH1("s2",xs2,0.2)
        th.FillH2("s1s2",xs1,xs2)

    print "histograms ={0}".format(th)
    th.DrawList(["s1","s2"],xd=1,yd=2)
    th.DrawAll()

    