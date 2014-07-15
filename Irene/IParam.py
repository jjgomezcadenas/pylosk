"""
Define running parameters
"""
from KParam import *

import os

generator ='Irene'  # can also be Toy, Paolina, etc
#inputFile="MagBox_Xe_10atm_00tesla.e2447.0.next"
inputFile="MagBox_Xe_10atm_00tesla.e2447.next"  # beware! set pressure right!
#inputFile="MagBox_Xe_10atm_00tesla.Xe136_bb0nu.next"
inputDir=os.environ['DATA']
outputDir='.'
pathToFile = inputDir+'/'+inputFile
debug = 0
gsEvtTree = "EVENT"
gsEvtBranch = "EventBranch"
eps = 0.001

Pr = 10.1325      #pressure in bar
sample =5.
Chi2Lim=4.

sigma_x=5. #measurement resolution in mm
sigma_y=5.


sigma_ux=0.3
sigma_uy=0.3

n_events=10000


DrawRoot = True
DrawMPL = True



	
