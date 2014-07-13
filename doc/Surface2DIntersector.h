#ifndef Surface2DIntersector_h
#define Surface2DIntersector_h 1

#include <recpack/ISurfaceIntersector.h>
#include <recpack/IModelEquation.h>
#include <recpack/LengthPropagator.h>

using namespace Recpack;

namespace Recpack{

  //! intersector between the helix state and a plane
  class  Surface2DIntersector: public ISurfaceIntersector  {
  public:
 
    //! default constructor
    Surface2DIntersector();
    
    //! default destructor
    virtual ~Surface2DIntersector(){}
    
    //! Matrix (ds/dx0 s=path length, x0= vector parameters at origin)
    const EMatrix& ds_dx0(const State& state, const Surface& surf, double length);

    
  protected:

    //! return equation
    IModelEquation& equation() 
      {return tool<IModelEquation>("equation");}

    //! return propagator
    LengthPropagator& propagator() 
      {return tool<LengthPropagator>("propagator");}


  };
}
#endif
  
