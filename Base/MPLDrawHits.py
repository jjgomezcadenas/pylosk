"""
MatPlotLib Implementation of Base class DrawHits
"""

from KDrawHits import KDrawHits
import matplotlib.pyplot as plt
import os
import sys
import numpy as np

from mpl_toolkits.mplot3d import Axes3D

class MPLDrawHits(KDrawHits):
    """
    Define an interface to draw hits
    """

    def __init__(self,binLimits,measurements,hits):
        """
        binLimits =([xl,yl,zl],[xu,yu,zu]]
        measurements = [x,y,z,e] or numpy array
        """

        KDrawHits.__init__(self,binLimits,measurements,hits)
        
        self.xl = self.lbin[0]
        self.yl = self.lbin[1]
        self.zl = self.lbin[2]
        self.xu = self.ubin[0]
        self.yu = self.ubin[1]
        self.zu = self.ubin[2]

        self.Meas = measurements
        self.Hits = hits

        self.meas_x0=[]
        self.meas_y0=[]
        self.meas_z0=[]
        self.meas_e0=[]

        self.hit_x0=[]
        self.hit_y0=[]
        self.hit_z0=[]
        self.hit_e0=[]
    
        for meas in self.Meas:
            self.meas_x0.append(meas[0])
            self.meas_y0.append(meas[1])
            self.meas_z0.append(meas[2])
            self.meas_e0.append(meas[3])

        for hit in self.Hits:
            self.hit_x0.append(hit[0])
            self.hit_y0.append(hit[1])
            self.hit_z0.append(hit[2])
            self.hit_e0.append(hit[3])
        


    def DrawMeasurements(self, draw='2D'):
        """
        Draws Measurements: the field draw selects 2D or 3D.
        """
        if draw == '2D':

            fig2D = plt.figure();
            #fig2D.set_figheight(15.0);
            #fig2D.set_figwidth(10.0);

            # Create the z-x projection.
            ax3 = fig2D.add_subplot(211)
            ax3.set_xlabel("z (cm)")
            ax3.set_ylabel("x (cm)")  
            ax3.set_xlim(self.zl, self.zu)
            ax3.set_ylim(self.xl, self.xu)  
    
            # Create the y-z projection.
            ax4 = fig2D.add_subplot(212)
            ax4.set_xlabel("z (cm)");
            ax4.set_ylabel("y (cm)");
            ax4.set_xlim(self.zl, self.zu)
            ax4.set_ylim(self.yl, self.yu)
            #self.ax2.plot(self.meas_x0,self.meas_y0,'o')
            ax3.plot(self.meas_z0,self.meas_x0,color='red')
            ax4.plot(self.meas_z0,self.meas_y0,color='red')
            plt.show()
            
        else:
            fig3D = plt.figure();
            #fig3D.set_figheight(15.0);
            #self.fig3D.set_figwidth(10.0);
            ax1 = fig3D.add_subplot(111, projection='3d');
            ax1.set_xlabel("x (cm)");
            ax1.set_ylabel("y (cm)");
            ax1.set_zlabel("z (cm)");
            ax1.set_xlim(self.xl, self.xu)
            ax1.set_ylim(self.yl, self.yu)
            ax1.set_zlim(self.zl, self.zu)

            ax1.plot(self.meas_x0,self.meas_y0,self.meas_z0,color='red')
            plt.show()

    
    

    def DrawHits(self, draw='2D'):
        """
        Draws hits: the field draw selects 2D or 3D.
        """

        if draw == '2D':
            fig2D = plt.figure();
            #fig2D.set_figheight(15.0);
            #fig2D.set_figwidth(10.0);

            # Create the x-z projection.
            ax3 = fig2D.add_subplot(211)
            ax3.set_xlabel("z (cm)")
            ax3.set_ylabel("x (cm)")  
            ax3.set_xlim(self.zl, self.zu)
            ax3.set_ylim(self.xl, self.xu)  

            # Create the y-z projection.
            ax4 = fig2D.add_subplot(212)
            ax4.set_xlabel("z (cm)");
            ax4.set_ylabel("y (cm)");
            ax4.set_xlim(self.zl, self.zu)
            ax4.set_ylim(self.yl, self.yu)
            
            ax3.plot(self.hit_z0,self.hit_x0,color='black');
            ax4.plot(self.hit_z0,self.hit_y0,color='black');
            plt.show()
        else:
            fig3D = plt.figure();
            #fig3D.set_figheight(15.0);
            #self.fig3D.set_figwidth(10.0);
            ax1 = fig3D.add_subplot(111, projection='3d');
            ax1.set_xlabel("x (cm)");
            ax1.set_ylabel("y (cm)");
            ax1.set_zlabel("z (cm)");
            ax1.set_xlim(self.xl, self.xu)
            ax1.set_ylim(self.yl, self.yu)
            ax1.set_zlim(self.zl, self.zu)
            ax1.plot(self.hit_x0,self.hit_y0,self.hit_z0,color='black')
            plt.show()

    
    def DrawAll(self, draw='2D'):
        """
        Draws measurements & hits: the field draw selects 2D or 3D.
        """
        if draw == '2D':
            fig2D = plt.figure();
            # Create the x-z projection.
            ax3 = fig2D.add_subplot(211)
            ax3.set_xlabel("z (cm)")
            ax3.set_ylabel("x (cm)")  
            ax3.set_xlim(self.zl, self.zu)
            ax3.set_ylim(self.xl, self.xu)  

            # Create the y-z projection.
            ax4 = fig2D.add_subplot(212)
            ax4.set_xlabel("z (cm)");
            ax4.set_ylabel("z (cm)");
            ax4.set_xlim(self.zl, self.zu)
            ax4.set_ylim(self.yl, self.yu)

            #ax2.plot(self.meas_x0,self.meas_y0,'o');
            ax3.plot(self.meas_z0,self.meas_x0,color='red');
            ax4.plot(self.meas_z0,self.meas_y0,color='red');
            #ax2.plot(self.hit_x0,self.hit_y0,color='black');
            ax3.plot(self.hit_z0,self.hit_x0,color='black');
            ax4.plot(self.hit_z0,self.hit_y0,color='black');
            plt.show()
        else:
            fig3D = plt.figure();
            ax1 = fig3D.add_subplot(111, projection='3d');
            ax1.set_xlabel("x (cm)");
            ax1.set_ylabel("y (cm)");
            ax1.set_zlabel("z (cm)");
            ax1.set_xlim(self.xl, self.xu)
            ax1.set_ylim(self.yl, self.yu)
            ax1.set_zlim(self.zl, self.zu)
            ax1.plot(self.meas_x0,self.meas_y0,self.meas_z0,color='red')
            ax1.plot(self.hit_x0,self.hit_y0,self.hit_z0,color='black')

    

    
    