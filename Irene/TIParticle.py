"""
Implementation of KParticle and Iparticle with concrete Root histograms
JJ, Spring, 2014
"""

from IParticle import IParticle
from ROOT import TH2F,TH3F
from ROOT import TCanvas, TPad, TPaveLabel, TPaveText, TColor
from ROOT import gROOT, gStyle
gStyle.SetOptStat(0);
gStyle.SetPalette(1);
gStyle.SetCanvasColor(33);
gStyle.SetFrameFillColor(18);

class TIParticle(IParticle):
    """
    Implements the interface to KMCParticle with an Irene particle
    """

    def __init__(self,ipart,itrk,smearVector,sample=5, Pressure=1):
        """
        Init the IParticle
        """

        IParticle.__init__(self,ipart,itrk,smearVector,sample,Pressure)
        
        self.__BookTrueHistograms()
        self.__FillTrueHistograms()
        self.__BookSmearHistograms()
        self.__FillSmearHistograms()

    def DrawHits(self, draw='2D',view='True'):
        """
        Draws hits: the field draw selects 2D or 3D.
        The field view can be equal to True for True hits, Smear for Smeared hits or both
        """
        
        if view == 'True':
            self.__DrawTrueHits(draw)
        elif view == 'Smear':
            self.__DrawSmearHits(draw)
        else:
            self.__DrawBothHits(draw)
    
    def __DrawTrueHits(self, draw='2D'):
        """
        Draw True hits: Drawing can be 2D or 3D
        """
        if draw == '2D':
            return self.__TDrawHits2D()
        else:
            return self.__TDrawHits3D()

    
    def __DrawSmearHits(self,draw='2D'):
        """
        Draw Smeared hits: Drawing can be 2D or 3D
        """

        if draw == '2D':
            return self.__TDrawSHits2D()
        else:
            return self.__TDrawSHits3D()


    def __DrawBothHits(self,draw='2D'):
        """
        Draw Smeared hits: Drawing can be 2D or 3D
        """

        if draw == '2D':
            return self.__TDrawBHits2D()
        else:
            return self.__TDrawBHits3D()
     

    def __TDrawHits2D(self):
        """
        Root implementation of 2D  
        """
        c1 = TCanvas( 'c1', 'True 2D', 200, 10, 600, 800 )
        c1.Divide(1,2)
        c1.cd(1)
        self.txz.Draw("colz")
        c1.cd(2)
        self.tyz.Draw("colz")

        s=raw_input("return to continue")

    def __TDrawHits3D(self):
        """
        Root implementation of 3D
        """
        c1 = TCanvas( 'c1', 'True 3D', 200, 10, 600, 800 )
        self.txyz.Draw("box")
        
        s=raw_input("return to continue")

    def __TDrawSHits2D(self):
        """
        Root implementation of 2D  
        """
        c1 = TCanvas( 'c1', 'Smear 2D', 200, 10, 600, 800 )
        c1.Divide(1,2)
        c1.cd(1)
        self.sxz.Draw("colz")
        c1.cd(2)
        self.syz.Draw("colz")

        s=raw_input("return to continue")

    def __TDrawSHits3D(self):
        """
        Root implementation of 3D
        """
        c1 = TCanvas( 'c1', 'Smear 3D', 200, 10, 600, 800 )
        self.sxyz.Draw("box")
        
        s=raw_input("return to continue")

    def __TDrawBHits2D(self):
        """
        Root implementation of 2D  
        """

        self.txz.SetMarkerColor(TColor.kRed)
        self.sxz.SetMarkerColor(TColor.kBlue)
        self.tyz.SetMarkerColor(TColor.kRed)
        self.syz.SetMarkerColor(TColor.kBlue)
        c1 = TCanvas( 'c1', 'Both 2D', 200, 10, 600, 800 )
        c1.Divide(1,2)
        c1.cd(1)
        self.txz.Draw("box")
        self.sxz.Draw("box same")
        c1.cd(2)
        self.tyz.Draw("box")
        self.syz.Draw("box same")

        s=raw_input("return to continue")

    def __TDrawBHits3D(self):
        """
        Root implementation of 3D
        """
        c1 = TCanvas( 'c1', 'Both 3D', 200, 10, 600, 800 )
        self.txyz.Draw("box")
        self.sxyz.Draw("box same")
        
        s=raw_input("return to continue")

    def __BookTrueHistograms(self):
        """
        Book the true histograms used for drawing
        """

        nhit = len(self.TrueHits())
        scale = self.Pr*0.5
        lbin = -700/scale
        ubin = 700/scale

        self.txz = TH2F("txz", "True: x vs z", nhit, lbin, ubin, nhit, lbin, ubin)
        self.tyz = TH2F("tyz", "True: y vs z", nhit, lbin, ubin, nhit, lbin, ubin)
        self.txyz = TH3F("tyxz", "True: x vs y vs z", nhit/50, lbin, ubin, 
                                                nhit/50, lbin, ubin,
                                                nhit/50, lbin, ubin)
        return 

    def __FillTrueHistograms(self):
        """
        Fill the true histograms used for drawing
        """

        for hit in self.TrueHits():
             self.txz.Fill(hit[2], hit[0], hit[3])
             self.tyz.Fill(hit[2], hit[1], hit[3])
             self.txyz.Fill(hit[0], hit[1], hit[2], hit[3])

        return 

    def __BookSmearHistograms(self):
        """
        Book the smear histograms used for drawing
        """

        nhit = len(self.SmearedHits())
        scale = self.Pr*0.5
        lbin = -700/scale
        ubin = 700/scale
        

        self.sxz = TH2F("sxz", "Smeared: x vs z", nhit, lbin, ubin, nhit, lbin, ubin)
        self.syz = TH2F("syz", "Smeared: y vs z", nhit, lbin, ubin, nhit, lbin, ubin)
        self.sxyz = TH3F("syxz", "Smeared: x vs y vs z", nhit/10, lbin, ubin, 
                                                nhit/10, lbin, ubin,
                                                nhit/10, lbin, ubin)
        return 

    def __FillSmearHistograms(self):
        """
        Fill the Smear histograms used for drawing
        """

        for hit in self.SmearedHits():
             self.sxz.Fill(hit[2], hit[0], hit[3])
             self.syz.Fill(hit[2], hit[1], hit[3])
             self.sxyz.Fill(hit[0], hit[1], hit[2], hit[3])

        return 

    