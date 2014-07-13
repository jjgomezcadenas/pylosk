"""
Module KFBase defines some basic objects for KF
"""
import numpy as np
import sys
eps=1e-10




class KFMatrix(object):
    """
    A KF matrix contains a numpy matrix
    """

    def __init__(self,obj):
            
        """
        M is a numpy array containing the coordinates
        """

        if isinstance(obj,list) ==True:
            self.M = np.matrix(obj)
        elif isinstance(obj,KFMatrix) ==True or isinstance(obj,KFVector) ==True:
            i,j = np.shape(obj.M)
            m=np.zeros(i*j)
            self.M=m.reshape(i,j)
            np.copyto(self.M,obj.M)
        else:
            print """input matrix must be a list with the following format
                    [[a11,a12,a13...],[a21,a22,a23...],[a31,a32,a33...]]
                    or another KFMatrix
                  """
            sys.exit(-1)


    def __eq__(self,other):
        m= abs(self.M-other.M)

        for i in range(0,np.shape(m)[0]):
            for j in range(0,np.shape(m)[1]):
                if m[i,j] > eps:
                    return False
        
        return True


    def __setitem__(self, (i, j), value): 
        """
        Allows m[i,j]=value 
        """ 
        
        self.M[i,j] = value 

    def __getitem__(self, (i, j)): 
        """
        Allows return m[i,j]   
        """ 

        return self.M[i,j]

    def __add__(self, other): 

        """
        M+x adds x to all elements of M
        M+N adds elements of M and N
        """

        
        if isinstance(other,KFMatrix) == True or isinstance(other,KFVector) == True:
            m = KFMatrix(self)
            m.M += other.M
            return m

        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m = KFMatrix(self)
            m.M += other
            return m
        else:
            print "cannot add a KFMatrix to type ",type(other)
            sys.exit(-1) 
        
    def __radd__(self, other): 

        """
        x+M
        """
        return self.__add__(other)

    def __sub__(self, other): 

        """
        M+x subtracts x to all elements of M
        M+N subtracts elements of M and N
        """

        
        if isinstance(other,KFMatrix) == True or isinstance(other,KFVector) == True:
            m = KFMatrix(self)
            m.M -= other.M
            return m

        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m = KFMatrix(self)
            m.M -= other
            return m
        else:
            print "cannot add a KFMatrix to type ",type(other)
            sys.exit(-1) 
        
    def __rsub__(self, other): 

        """
        x-M
        """
        return self.__sub__(other)

    def __mul__(self, other): 

        """
        M*Y matrix multiplication if both are matrices
        M*e multiplication of element if e a number
        """
        
        if isinstance(other,KFMatrix) == True:
            m = KFMatrix(self)
            m.M = np.dot(self.M,other.M)
            return m
        elif isinstance(other,KFVector) == True:
            m = KFVector(other)
            m.M = np.dot(self.M,other.M)
            return m
        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m = KFMatrix(self)
            m.M *= other
            return m
        else:
            print "cannot multiply a KFMatrix to type ",type(other)
            sys.exit(-1) 


    def __rmul__(self, other): 

        """
        M*Y matrix multiplication if both are matrices
        M*e multiplication of element if e a number
        """

        m = KFMatrix(self)
        if isinstance(other,KFMatrix) == True or isinstance(other,KFVector) == True:
            m.M = np.dot(other.M,self.M)

        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m.M *= other
            
        else:
            print "cannot multiply a KFMatrix to type ",type(other)
            sys.exit(-1) 

        return m

    def Shape(self):
        return np.shape(self.M)

    def Transpose(self):
        m = KFMatrix(self)
        m.M=self.M.T
        return m

    def Determinant(self):
        return np.linalg.det(self.M)

    def Inverse(self):
        m = KFMatrix(self)
        try:
            m.M=np.linalg.inv(self.M)
        except np.linalg.linalg.LinAlgError:
            print "WARNING: singular matrix, cannot properly invert..."
        return m

    def __str__(self):
        return self.M.__str__()

    def __repr__(self):
        return self.M.__repr__()


class KFVector(object):
    """
    A KF vector contains a numpy array
    """

    def __init__(self,obj):
        """
        Takes as an argument a list or another KFVector.
        Represents the vector as a column 
        """
        if isinstance(obj,list) ==True:
            self.M = np.transpose(np.matrix(obj))
            row,col = np.shape(self.M)
            if not col == 1:
                print "A KFVector takes as input a list of number"
                print "for example: v = KFVector([1,2,3])"
                sys.exit(-1)

        elif isinstance(obj,KFVector) ==True or isinstance(obj,KFMatrix) ==True:
            i,j = np.shape(obj.M)
            if j != 1:
                print "A KFVector must be a column vector"
                print "Imput was:", obj.__str__()
                sys.exit(-1)

            m=np.zeros(i*j)
            self.M=m.reshape(i,j)
            np.copyto(self.M,obj.M)
            self.M=self.M.reshape(i*j,1)
        
        else:
            print """input  must be a list v = KFVector([1,2,3])
                     or another KFVector
                  """
            sys.exit(-1)



    def __eq__(self,other):
        m= abs(self.M-other.M)

        for i in range(0,np.shape(m)[0]):
            for j in range(0,np.shape(m)[1]):
                if m[i,j] > eps:
                    return False
        
        return True

    def __setitem__(self, i, value): 
        """
        Allows m[i]=value 
        """ 
        
        self.M[i,0] = value 

    def __getitem__(self, i): 
        """
        Allows return m[i]   
        """ 

        return self.M[i,0]


    def __add__(self, other): 

        """
        M+x adds x to all elements of M
        M+N adds elements of M and N
        """

        m = KFVector(self)

        if isinstance(other,KFVector) ==True or isinstance(other,KFMatrix) ==True:
            row,col = np.shape(other.M)
            if not col == 1:
                print "A KFVector must be a column vector"
                print "Imput was:", other.__str__()
                sys.exit(-1)

            m.M += other.M
        
        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m.M += other 
        else:
            print "cannot add the KFVector to type ",type(other)
            sys.exit(-1)
            

        return m

    def __sub__(self, other): 

        """
        M+x subtracts x to all elements of M
        M+N subtracts elements of M and N
        """

        m = KFVector(self)

        if isinstance(other,KFVector) ==True or isinstance(other,KFMatrix) ==True:
            row,col = np.shape(other.M)
            if not col == 1:
                print "A KFVector must be a column vector"
                print "Imput was:", other.__str__()
                sys.exit(-1)

            m.M -= other.M
        
        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m.M -= other 
        else:
            print "cannot add the KFVector to type ",type(other)
            sys.exit(-1)
            

        return m


    def __radd__(self, other): 

        """
        x+M
        """
        return self.__add__(other)

    def __rsub__(self, other): 

        """
        x+M
        """
        return self.__sub__(other)

    def __mul__(self, other): 

        """
        product of two KFVectors ---> scalar product (dot)
        Product by a matrix H*a ---> a KFVector
        Product by a matrix a*H ---> a KFMatrix
        """

        if isinstance(other,KFMatrix) == True:
            m = KFMatrix(self)
            m.M = np.dot(self.M,other.M)
            return m
        elif isinstance(other,KFVector) == True:
            return np.dot(np.transpose(other.M),self.M)
        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m = KFVector(self)
            m.M*=other
            return m
        else:
            print "cannot multiply a KFVector to type ",type(other)
            sys.exit(-1) 

    def __rmul__(self, other): 

        """
        product of two KFVectors ---> scalar product (dot)
        Product by a matrix H*a ---> a KFVector
        Product by a matrix a*H ---> a KFMatrix
        """

        if isinstance(other,KFMatrix) == True:
            m = KFVector(self)
            m.M = other.M * self.M 
            return m
        elif isinstance(other,KFVector) == True:
            return np.sum(np.transpose(other.M)*self.M)

        elif isinstance(other,int) ==True or isinstance(other,float) == True: 
            m = KFVector(self)
            m.M *= other
            return m
        else:
            print "cannot multiply a KFVector to type ",type(other)
            sys.exit(-1)

    def Shape(self):
        return np.shape(self.M)

    def Transpose(self):
        m = KFMatrix(self)
        m.M=self.M.T
        return m

    def Length(self):
        return len(self.M)

    
    def Mod(self):
        return np.linalg.norm(self.M)

    def Norm(self):
        return np.linalg.norm(self.M)

    def Sum(self):
        return np.sum(self.M)

    def __str__(self):
        return self.M.__str__()

    def __repr__(self):
        return self.M.__repr__()


class KFHVector(object):
    """

    A KHVector is built from a column vector of dimension (1,p) 
    and a matrix of dimension (p,p).
    It is the base clase to represent 
    A measurement characterized by a hit (for example (x,y)) and
    a Covariance matrix. 
    
    A state, characterized by a trajectory (x,y,ux,uy) and its Covariance Matrix
    
    """
    
    def __init__(self,vect, cov):
        """
        vect: a KFVector of measurements (eg, x,y) or trajectory params (x,y,ux,uy) 
        cov a KFMatrix of Covariance (
        
        """

        if isinstance(vect, KFVector) != True and isinstance(vect, KFMatrix) != True:
            print "the vector must be a KFVector (or a KFMatrix)"
            sys.exit(-1)

        if isinstance(cov, KFMatrix) != True:
            print "the Covariance must be a KFMatrix"
            sys.exit(-1)

      
        self.V = vect
        self.Cov = cov


    def __str__(self):
        s ="<\nVector :"
        s+=self.V.__str__()+"\nCov :"+self.Cov.__str__()+"\n>"
        return s 
    def __repr__(self):
        return self.__str__()  

def testVector():
    v1 = KFVector([1.,2.,3.])  #notice: 1.,2.,3.] to get floats
    vv2 = KFVector([2.,3.,4.]) 
    v1pv2=KFVector([3.,5.,7.])
    v1mv2=KFVector([-1.,-1.,-1.])
    v2mv1=KFVector([1.,1.,1.])

    c1 = v1.Shape()==(3,1)

    print "test property: Shape()"
    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected v1.Shape() = (3,1)
        found ={0}
        """.format(v1.Shape())
        sys.exit(-1)

    c1 = v1.Length()==3

    print "test property: Length()"
    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected v1.Length() = 3
        found ={0}
        """.format(v1.Length())
        sys.exit(-1)

    c1 = abs(v1.Norm()-3.74165738677) < 1e-6

    print "test property: Norm()"
    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected v1.Norm() = 3.74165738677
        found ={0}
        """.format(v1.Norm())
        sys.exit(-1)

    c1 = v1[0]==1. and v1[1]==2. and v1[2]==3.

    print "test property: access to element"
    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected v1[0]==1. and v1[1]==2. and v1[2]==3.
        found ={0}
        """.format(v1)
        sys.exit(-1)

    v2 = KFVector(v1)
    c1 = v2.Shape()==(3,1) and v2.Length()==3 
    c1 = c1 and abs(v2.Norm()-3.74165738677) < 1e-6 
    c1 = c1 and v2[0]==1. and v2[1]==2. and v2[2]==3.

    print "test property: Create vector from vector"

    if c1 == True:
        print "passed"
    else:
        print """
        failed: 
        vector v2 = {0}
        shape ={1}
        length ={2}
        norm ={3}
        type ={4}
        """.format(v2,v2.Shape(),v2.Length(), v2.Norm(),type(v2))
        sys.exit(-1)

    print "test property: v1 = v2"

    c1 = v1 ==v2
    if c1 == True:
        print "passed"
    else:
        print """
        failed
    
        v1 ={0}
        type(v1) ={1}
        v2 = {2}
        type(v2) ={3} 
        """.format(v1,type(v1),v2,type(v2))
        sys.exit(-1)

    print "test property: write access to element"

    for i in range(0,v2.Length()):
        v2[i]+=1.

    c1 = v2[0]==2. and v2[1]==3. and v2[2]==4.

    
    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected v2[0]==2. and v2[1]==3. and v2[2]==4.
        found ={0}
        """.format(v2)
        sys.exit(-1)

    print "test property: v1 != v2"

    c1 = v1 !=v2
    if c1 == True:
        print "passed"
    else:
        print """
        failed
    
        v1 ={0}
        type(v1) ={1}
        v2 = {2}
        type(v2) ={3} 
        """.format(v1,type(v1),v2,type(v2))
        sys.exit(-1)

    print "test property: v1 * v1"

    c1 = v1*v1 == 14

    if c1 == True:
        print "passed"
    else:
        print """
        v1*v1
        ------
        type ={0}
        v1*v1 = {1}
        """.format(type(v1*v1),v1*v1)
        sys.exit(-1)

    print "test property: v1 * v2"

    c1 = v1*v2 == v2*v1 == 20

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        v1 ={0}
        type(v1) ={1}
        v2 = {2}
        type(v2) ={3} 
        """.format(v1,type(v1),v2,type(v2))
        sys.exit(-1)

    print "test property: v1 + v2"

    v1pv2=KFVector([3.,5.,7.])

    c1 = v1+v2 == v2+v1 == v1pv2

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        v1 ={0}
        type(v1) ={1}
        v2 = {2}
        type(v2) ={3} 
        """.format(v1,type(v1),v2,type(v2))
        sys.exit(-1)

    print "test property: v1 - v2"

    c1 = v1-v2 == v1mv2

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        v1-v2 ={0}
        v1mv2={1} 
        """.format(v1-v2,v1mv2)
        sys.exit(-1)

    print "test property: v2 - v1"

    c1 = v2-v1 == v2mv1

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        v2-v1 ={0}
        v2mv1={1} 
        """.format(v2-v1,v2mv1)
        sys.exit(-1)


    print "test property: v1 * number"

    v1p10=KFVector([10.,20.,30.])

    c1 = v1*10 == 10*v1 == v1p10

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        v1 ={0}
        type(v1) ={1}
        v1*10 = {2}
        type(v1*10) ={3} 
        """.format(v1,type(v1),v1*10,type(v1*10))
        sys.exit(-1)

    print "test property: v1 + number"

    v1p10=KFVector([11.,12.,13.])

    c1 = v1+10 == 10+v1 == v1p10

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        v1 ={0}
        type(v1) ={1}
        v1+10 = {2}
        type(v1+10) ={3} 
        """.format(v1,type(v1),v1+10,type(v1+10))
        sys.exit(-1)


def exampleVector():

    v1 = KFVector([1.,2.,3.])  #notice: 1.,2.,3.] to get floats

    print """
    Create vector from list
    -----------------------

    vector v1 = {0}
    shape ={1}
    length ={2}
    norm ={3}
    type ={4}
    """.format(v1,v1.Shape(),v1.Length(), v1.Norm(),type(v1))

    print """
             v1[0] = {0} 
             v1[1] = {1}  
             v1[2] = {2}     
               """.format(v1[0],v1[1],v1[2])

    v2 = KFVector(v1)    

    print """
    Create vector from vector
    -------------------------

    vector v2 = {0}
    shape ={1}
    length ={2}
    norm ={3}
    type ={4}
    """.format(v2,v2.Shape(),v2.Length(), v2.Norm(),type(v2))

    print """
    Is v1 = v2?
    -----------
    """

    if v1 == v2:
        print "True"
    else:
        print "False"

    print """
    Change v2 elements
    -------------------
    """
    for i in range(0,v2.Length()):
        v2[i]+=1.

    print """
    v1 ={0}
    v2 ={1}
    """.format(v1,v2)

    print """
    Is v1 = v2?
    -----------
    """

    if v1 == v2:
        print "True"
    else:
        print "False"


    print """
    v1*v1
    ------
    type ={0}
    v1*v1 = {1}

    """.format(type(v1*v1),v1*v1)

    print "is numpy float a float?",isinstance(v1*v1,float)

    norm = v1.Norm()
    vnorm = v1*(1./norm)

    print """
    Norm
    -------
    norm ={0}
    1./norm = {1}
    
    """.format(norm,1./norm)

    print """
    vnorm
    -------
    vnorm = {0}
    shape ={1}
    length ={2}
    norm ={3}
    Sum(vnorm) ={4}
    """.format(vnorm,vnorm.Shape(),vnorm.Length(), vnorm.Norm(), vnorm.Sum())

    print """
    v1 & vnorm
    ---------
    v1 ={0}
    type(v1) ={1}
    vnorm = {2}
    type(vnorm) ={3} 
    """.format(v1,type(v1),vnorm,type(vnorm))

    print """
    v1*vnorm
    -----
    type ={0}
    v1*v2 = {1}

    vnorm*v1
    -----
    type ={2}
    v2*v1 = {3}

    """.format(type(v1*vnorm),v1*vnorm,type(vnorm*v1),vnorm*v1)

    print """
    v1 & v2
    ---------
    v1 ={0}
    type(v1) ={1}
    v2 = {2}
    type(v2) ={3} 
    """.format(v1,type(v1),v2,type(v2))

    print """
    v1*v2
    -----
    type ={0}
    v1*v2 = {1}

    v2*v1
    -----
    type ={2}
    v2*v1 = {3}

    """.format(type(v1*v2),v1*v2,type(v2*v1),v2*v1)


    print """
    v1 & v2
    ---------
    v1 ={0}
    v2 = {1} 
    """.format(v1,v2)

    print """
    v1+v2
    -------
    type ={0}
    v1+v2 = {1}

    v2+v1
    ------
    type ={2}
    v2+v1 = {3}

    """.format(type(v1+v2),v1+v2,type(v2+v1),v2+v1)

    print """
    v1 & v2
    ---------
    v1 ={0}
    v2 = {1} 
    """.format(v1,v2)

    print """
    v1*number
    ---------
    v1 ={0}
    v1*10 = {1}
    10*v1= {2} 
    """.format(v1,v1*10,10*v1)

    print """
    v1+number
    v1+10 = {0}
    10+v1= {1} 
    """.format(v1+10,10+v1)


def exampleMatrix():
             
    m =KFMatrix([[1.,2.,3.],[0.,1.,4.],[5.,6.,0.]])
    v1 = KFVector([1.,2.,3.])  #notice: 1.,2.,3.] to get floats

    print """
    matrix m = {0}
    shape ={1}
    Transpose ={2}
    Determinant ={3}
    Inverse ={4}
    type = {5}

    """.format(m,m.Shape(),m.Transpose(), m.Determinant(), m.Inverse(),type(m))

                

    print """
             m[0,0] = %7.1f  m[0,1] = %7.1f  m[0,2] = %7.1f 
             m[1,0] = %7.1f  m[1,1] = %7.1f  m[1,2] = %7.1f  
             m[2,0] = %7.1f  m[2,1] = %7.1f  m[2,2] = %7.1f   
               """%(
             m[0,0],  m[0,1] ,  m[0,2] ,   
             m[1,0],  m[1,1] ,  m[1,2] ,   
             m[2,0],  m[2,1] ,  m[2,2] ,   
                )
    print """
    m*v1
    type ={0}
    m*v1 = {1}

    """.format(type(m*v1),m*v1)

    print """
    matrix m*m^-1 
    m ={0}
    m^-1 ={1}
    m*m^-1 ={2}
    m^-1*m ={3}

    """.format(m,m.Inverse(),m*m.Inverse(),m.Inverse()*m)

    print """
    m*number
    m*10 = {0}
    10*m= {1} 
    """.format(m*10,10*m)

    print """
    m+number
    m+10 = {0}
    10+m= {1} 
    """.format(m+10,10+m)


def testMatrix():
    v1 = KFVector([1.,2.,3.])  #notice: 1.,2.,3.] to get floats

    m =KFMatrix([[1.,2.,3.],[0.,1.,4.],[5.,6.,0.]])
    mp10 =KFMatrix([[10.,20.,30.],[0.,10.,40.],[50.,60.,0.]])
    ma10 =KFMatrix([[11.,12.,13.],[10.,11.,14.],[15.,16.,10.]])
    mm =KFMatrix([[16.,22.,11.],[20.,25.,4.],[5.,16.,39.]])
    mT =KFMatrix([[1.,0.,5.],[2.,1.,6.],[3.,4.,0.]])
    mI=KFMatrix([[-24.,  18.,   5.],[20., -15.,  -4.],[-5.,   4.,   1.]])
    mpmT=KFMatrix([[ 2.,  2.,  8.],
        [ 2.,  2., 10.],
        [ 8., 10.,  0.]])

    I =KFMatrix([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])

    print """
    matrix m = {0}
    shape ={1}
    Transpose ={2}
    Determinant ={3}
    Inverse ={4}
    type = {5}

    """.format(m,m.Shape(),m.Transpose(), m.Determinant(), m.Inverse(),type(m))

                
    print "test property: Shape()"
    c1 = m.Shape()==(3,3)

    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected v1.Shape() = (3,3)
        found ={0}
        """.format(m.Shape())
        sys.exit(-1)


    c1 = m[0,0]==1. and m[1,1]==1. and m[2,2]==0.

    print "test property: access to element"
    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected m[0,0]==1. and m[1,1]==1. and m[2,2]==0.
        found ={0}
        """.format(m)
        sys.exit(-1)

    m2 = KFMatrix(m)
    c1 = m2.Shape()==(3,3) 
    c1 = c1 and m2[0,0]==1. and m2[1,1]==1. and m2[2,2]==0.

    print "test property: Create matrix from matrix"

    if c1 == True:
        print "passed"
    else:
        print """
        failed: 
        matrix m2 = {0}
        shape ={1}
        type ={2}
        """.format(m2,m2.Shape(),type(m2))
        sys.exit(-1)

    print "test property: m = m2"

    c1 = m ==m2
    if c1 == True:
        print "passed"
    else:
        print """
        failed
    
        m ={0}
        type(m) ={1}
        m2 = {2}
        type(m2) ={3} 
        """.format(m,type(m),m2,type(m2))
        sys.exit(-1)

    print "test property: write access to element"

    for i in range(0,3):
        m2[i,i]+=1.

    c1 = m2[0,0]==2. and m2[1,1]==2. and m2[2,2]==1.

    
    if c1 == True:
        print "passed"
    else:
        print """
        failed: expected m2[0,0]==2. and m2[1,1]==2. and m2[2,2]==1.
        found ={0}
        """.format(m2)
        sys.exit(-1)

    print "test property: m != m2"

    c1 = m !=m2
    if c1 == True:
        print "passed"
    else:
        print """
        failed
    
        m ={0}
        type(m) ={1}
        m2 = {2}
        type(m2) ={3} 
        """.format(m,type(m),m2,type(m2))
        sys.exit(-1)

    print "test property: m * m"

    c1 = m*m == mm

    if c1 == True:
        print "passed"
    else:
        print """
        m*m
        ------
        m*m ={0}
        mm = {1}
        """.format(m*m,mm)

    print "test property: Transpose"

    c1 = m.Transpose() == mT 

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        
        m.Transpose() ={0}
        mT = {1}
        """.format(m.Transpose(),mT)

    print "test property: Inverse"

    c1 = m.Inverse() == mI 

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        
        m.Inverse() ={0}
        mI = {1}
        """.format(m.Inverse(),mI)

    print "test property: m*m^-1"

    c1 = m*m.Inverse() == m.Inverse()*m ==I 

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        
        m*m.Inverse() ={0}
        m.Inverse()*m ={1}
        I = {2}
        """.format(m*m.Inverse(),m.Inverse()*m,I)

    print "test property: m1 + m2"

    
    c1 = m+m.Transpose() == m.Transpose()+m == mpmT

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        m+m.Tranpose() ={0}
        m.Tranpose()+m ={1}
        mpmT= {2}
        """.format(m+m.Tranpose(),m.Tranpose()+m,mpmT)


    print "test property: v1 * number"

    c1 = m*10 == 10*m == mp10

    if c1 == True:
        print "passed"
    else:
        print """
        Failed
        ---------
        mp10={0}
        m*10 = {1}
        """.format(mp10,m*10)

    print "test property: v1 + number"

    c1 = m+10 == 10+m == ma10

    if c1 == True:
        print "passed"
    else:
        print """
        
        Failed
        ---------
        ma10 ={0}
        m+10 = {1}
        """.format(ma10,m+10)


def testKHVector():
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

    m = KFHVector(hit,cov)
    
    print "Measurement=",m

def testMatrixMult():
    
    print "\n\n** TESTING MATRIX MULTIPLICATION **\n\n"
    
    m1 = KFMatrix([[1.,2.,-3.,4.],[0.,1.,-4.,3.],[5.,-6.,0.,2.],[-2.,4.,1.,-5.]])
    m2 = KFMatrix([[2.,-5.,1.,-2.],[2.,1.,1.,2.],[4.,-3.,3.,2.],[1.,1.,6.,-3.]])
    v1 = KFVector([2.,3.,-6.,1.])
    v2 = KFVector([0.,2.,2.,4.])

    prod1_expect = KFMatrix([[-2.,10.,18.,-16.],[-11.,16.,7.,-15.],[0.,-29.,11.,-28.],[3.,6.,-25.,29.]]) # m1*m2
    prod2_expect = -2  # v1*v2
    prod3_expect = KFVector([30.,30.,-6.,-3.])
    prod4_expect = KFMatrix([-13.,12.,-7.,-13.])
    
    print "Computing m1*m2"
    prod1 = m1*m2
    c1 = (prod1 == prod1_expect)
    print prod1
    if(c1):
        print "\n\n-- PASSED --\n\n"
    else:
        print "\n\n-- FAILED --\n\n"
        
    print "Computing v1*v2"    
    prod2 = v1*v2
    c2 = (prod2 == prod2_expect)
    print prod2
    if(c2):
        print "\n\n-- PASSED --\n\n"
    else:
        print "\n\n-- FAILED --\n\n"

    print "Computing m1*v1"    
    prod3 = m1*v1
    print prod3
    c3 = (prod3 == prod3_expect)
    if(c3):
        print "\n\n-- PASSED --\n\n"
    else:
        print "\n\n-- FAILED --\n\n"
    
    print "Computing v1T*m2"
    prod4 = v1.Transpose()*m2
    c4 = (prod4 == prod4_expect)
    print prod4
    if(c4):
        print "\n\n-- PASSED --\n\n"
    else:
        print "\n\n-- FAILED --\n\n"    

if __name__ == '__main__':
    exampleVector()
    testVector()
    exampleMatrix()
    testMatrix()
    #testKHVector()
    testMatrixMult()
