"""
Root Implementation of Base class DrawHits
"""

from KDrawHits import KDrawHits
from ROOT import TH2F,TH3F
from ROOT import TCanvas, TPad, TPaveLabel, TPaveText, TColor
from ROOT import gROOT, gStyle
gStyle.SetOptStat(0);
gStyle.SetPalette(1);
gStyle.SetCanvasColor(33);
gStyle.SetFrameFillColor(18);

class TDrawHits(KDrawHits):
    """
    Define an interface to draw hits
    """

    def __init__(self,binLimits,measurements,hits):
        """
        binLimits =([xl,yl,zl],[xu,yu,zu]]
        measurements = [x,y,z,e] or numpy array
        """

        KDrawHits.__init__(self,binLimits,measurements,hits)
        
        
        self.__BookMeasHistograms()
        self.__FillMeasHistograms()
        self.__BookHitHistograms()
        


    def DrawMeasurements(self, draw='2D'):
        """
        Draws Measurements: the field draw selects 2D or 3D.
        """
        if draw == '2D':
            return self.__TDrawMeas2D()
        else:
            return self.__TDrawMeas3D()


    def DrawHits(self, draw='2D'):
        """
        Draws hits: the field draw selects 2D or 3D.
        """

        self.__FillHitHistograms()

        if draw == '2D':
            return self.__TDrawHits2D()
        else:
            return self.__TDrawHits3D()

    
    def DrawAll(self, draw='2D'):
        """
        Draws measurements & hits: the field draw selects 2D or 3D.
        """
        if draw == '2D':
            return self.__TDrawAll2D()
        else:
            return self.__TDrawAll3D()

    def __TDrawMeas2D(self):
        """
        Root implementation of 2D  
        """
        c1 = TCanvas( 'c1', 'Measurements', 200, 10, 600, 800 )
        c1.Divide(1,2)
        c1.cd(1)
        self.mxz.Draw("colz")
        c1.cd(2)
        self.myz.Draw("colz")

        s=raw_input("return to continue")

    def __TDrawMeas3D(self):
        """
        Root implementation of 3D
        """
        c1 = TCanvas( 'c1', 'Measurement', 200, 10, 600, 800 )
        self.mxyz.Draw("box")
        
        s=raw_input("return to continue")

    def __TDrawHits2D(self):
        """
        Root implementation of 2D  
        """
        c1 = TCanvas( 'c1', 'Hits', 200, 10, 600, 800 )
        c1.Divide(1,2)
        c1.cd(1)
        self.hxz.Draw("colz")
        c1.cd(2)
        self.hyz.Draw("colz")

        s=raw_input("return to continue")

    def __TDrawHits3D(self):
        """
        Root implementation of 3D
        """
        c1 = TCanvas( 'c1', 'Hits', 200, 10, 600, 800 )
        self.hxyz.Draw("box")
        
        s=raw_input("return to continue")

    def __TDrawAll2D(self):
        """
        Root implementation of 2D  
        """

        self.mxz.SetMarkerColor(TColor.kRed)
        self.hxz.SetMarkerColor(TColor.kBlue)
        self.myz.SetMarkerColor(TColor.kRed)
        self.hyz.SetMarkerColor(TColor.kBlue)
        self.mxz.SetMarkerStyle(6)
        #self.mxz.SetMarkerSize(1)
        self.hxz.SetMarkerStyle(7)
        #self.hxz.SetMarkerSize(2)
        self.myz.SetMarkerStyle(6)
        self.hyz.SetMarkerStyle(7)
        #self.hyz.SetMarkerSize(2)
        c1 = TCanvas( 'c1', 'Meas & Hits 2D', 200, 10, 600, 800 )
        c1.Divide(1,2)
        c1.cd(1)
        self.mxz.Draw("")
        self.hxz.Draw("same")
        c1.cd(2)
        self.myz.Draw("")
        self.hyz.Draw("same")

        s=raw_input("return to continue")

    def __TDrawAll3D(self):
        """
        Root implementation of 3D
        """
        c1 = TCanvas( 'c1', 'Both 3D', 200, 10, 600, 800 )
        self.mxyz.Draw("box")
        self.hxyz.Draw("box same")
        
        s=raw_input("return to continue")

    def __BookMeasHistograms(self):
        """
        Book the true histograms used for drawing
        """


        xl = self.lbin[0]
        yl = self.lbin[1]
        zl = self.lbin[2]
        xu = self.ubin[0]
        yu = self.ubin[1]
        zu = self.ubin[2]

        nhit = len(self.Meas)
        

        self.mxz = TH2F("mxz", "Meas: x vs z", nhit, zl, zu, nhit, xl, xu)
        self.myz = TH2F("myz", "Meas: y vs z", nhit, zl, zu, nhit, yl, yu)
        self.mxyz = TH3F("myxz", "Meas: x vs y vs z", nhit/50, xl, xu, 
                                                nhit/50, yl, yu,
                                                nhit/50, zl, zu)
        return 

    def __FillMeasHistograms(self):
        """
        Fill the true histograms used for drawing
        """

        for hit in self.Meas:
             self.mxz.Fill(hit[2], hit[0], hit[3])
             self.myz.Fill(hit[2], hit[1], hit[3])
             self.mxyz.Fill(hit[0], hit[1], hit[2], hit[3])

        return 

    def __BookHitHistograms(self):
        """
        Book the smear histograms used for drawing
        """

        xl = self.lbin[0]
        yl = self.lbin[1]
        zl = self.lbin[2]
        xu = self.ubin[0]
        yu = self.ubin[1]
        zu = self.ubin[2]

        nhit = len(self.Hits)
        

        self.hxz = TH2F("hxz", "Meas: x vs z", nhit, zl, zu, nhit, xl, xu)
        self.hyz = TH2F("hyz", "Meas: y vs z", nhit, zl, zu, nhit, yl, yu)
        self.hxyz = TH3F("hyxz", "Meas: x vs y vs z", nhit/10, xl, xu, 
                                                nhit/10, yl, yu,
                                                nhit/10, zl, zu)
        return 

    def __FillHitHistograms(self):
        """
        Fill the Smear histograms used for drawing
        """

        for hit in self.Hits:
             self.hxz.Fill(hit[2], hit[0], hit[3])
             self.hyz.Fill(hit[2], hit[1], hit[3])
             self.hxyz.Fill(hit[0], hit[1], hit[2], hit[3])

        return 

    
    