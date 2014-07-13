#ifndef HelixMsNoiser_h
#define HelixMsNoiser_h 1

#include <recpack/MsNoiser.h>
#include <recpack/IStateRepresentation.h>
#include <recpack/HelixEquation.h>
#include <recpack/Propagator.h>

using namespace Recpack;

namespace Recpack {

  //! multiple scattering error for a the straight line model
  /* NOTE: it will ask for the hypervector named "overp" in the state*/
  class HelixMsNoiser : public MsNoiser {

  public:

    //! default constructor with a radiation length X0 and position
    /* position it is the index where to locate the slopes matrix*/
    HelixMsNoiser(double x0, size_t ipos);

    //! default destructor
    virtual ~HelixMsNoiser() {};
    
    //! compute and return the Q matrix matrix 
    const EMatrix& QMatrix(const State& state, double length);
    
  protected:

    HelixEquation& equation()
      {return tool<HelixEquation>("equation");}

    IStateRepresentation& representation()
      {return tool<IStateRepresentation>("representation");}

    //! return state equation
    Propagator& propagator(){
      return tool<Propagator>("propagator");
    }

    MsNoiser& sline() 
      {return tool<MsNoiser>("noiser/sline/ms");}

  protected:

    //! position where to colocate the slopes
    size_t _i;

  };
  
}
#endif
