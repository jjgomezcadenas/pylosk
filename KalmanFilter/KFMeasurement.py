import sys
from KFBase import KFVector, KFMatrix, KFHVector

class KFMeasurement(KFHVector):
    """

    Represents a measurement characterized by a hit (for example (x,y)) and
    a Covariance matrix. 
    The hit is represented by a KFVector and the Covariance Matrix by a KFMatrix
    CovarianceMatrix=KFMatrix([[sigma_x**2,0],
                   [0,sigma_y**2]])
    """
    
    def __init__(self,hit, cov):
        """
        hit: a KFVector of measurements (eg, x,y) 
        cov a KFMatrix of Covariance (diagonal: sigma_x^2,sigma_y^2
        
        """

        KFHVector.__init__(self,hit,cov)
        

    def __str__(self):
        s ="<\nHit :"
        s+=self.V.__str__()+"\nCov :"+self.Cov.__str__()+"\n>"
        return s

    def __repr__(self):
        return self.__str__()   
    
if __name__ == '__main__':

    x0=5.
    y0=15.
    hit =KFVector([x0,y0])
    sigma_x=0.1
    sigma_y=0.2
    cov=KFMatrix([[sigma_x**2,0],
                   [0,sigma_y**2]])

    print "hit is a KFVector",isinstance(hit,KFVector)
    print "hit =",hit

    print "cov is a KFMatrix",isinstance(cov,KFMatrix)
    print "cov =",cov

    m = KFMeasurement(hit,cov)
    
    print "Measurement=",m
    
    
    
