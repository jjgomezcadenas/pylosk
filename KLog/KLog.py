"""
Define logging parameters
"""
import logging 
from KEnum import Enum

Debug=Enum("DEBUG", mute=1,quiet=2,info=3,verbose=4,draw=5)


#create Hanler
ch = logging.StreamHandler()
#ch.setLevel(logging.INFO)
#create formatter and add it to the handler
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger


def cond_pause(debug_level):
    """
    Pauses if the debug_level is smaller than the current debug level 
    """
    if debug_level >= Debug.verbose.value:
        s=raw_input("return to continue")
