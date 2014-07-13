#ifndef Helix_PlaneIntersector_h
#define Helix_PlaneIntersector_h 1

#include <recpack/Surface2DIntersector.h>
#include <recpack/HelixEquation.h>

using namespace Recpack;

namespace Recpack{

  //! intersector between the helix state and a plane
  class  Helix_PlaneIntersector: public Surface2DIntersector  {
  public:
 
    //! default constructor
    Helix_PlaneIntersector(); 
    
    //! default destructor
    virtual ~Helix_PlaneIntersector(){}
    
    //! compute the intersection and return length as argument
    bool length(const State& state, const Surface& surf, 
		double& s);

    //! Matrix (ds/dS s=path length, S= vector surface parameters -- position and direction --- )
    const EMatrix& ds_dS(const State& state, const Surface& surf, double length);
    
  protected:
    
    bool initialize(const State& state, const Surface& surf);
    
    bool caseBParallel(double& s);
    
    bool caseBPerpendicular(double& s);
    
    bool caseBAproximation(double& s);

    ISurfaceIntersector& numerical() 
      {return tool<ISurfaceIntersector>("intersector/numerical");}

    ISurfaceIntersector& helix_numerical() 
      {return tool<ISurfaceIntersector>("intersector/helix_num");}

    ISurfaceIntersector& sline() 
      {return tool<ISurfaceIntersector>("intersector/sline/plane");}

    
  protected:


    HelixEquation* _helix;

    EVector _rp;
    EVector _np;
    
    EVector _b;
    EVector _r0;
    EVector _u0;
    EVector _n0;
    EVector _m0;
    
    double _alpha;
    double _gamma;
    double _QQ;
    double _RR;
    double _charge;
    
  };
}
#endif
