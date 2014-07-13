#ifndef HelixPropagator_h
#define HelixPropagator_h 1

#include <recpack/LengthPropagator.h>
#include <recpack/HelixEquation.h>

using namespace Recpack;

namespace Recpack {

  //! propagato of the helix, a B field should be stored
  class HelixPropagator : public LengthPropagator {
  public:

    //! default constructor
    HelixPropagator() {_sdim=3;}

    //! default destructor
    virtual ~HelixPropagator() {}
    
  protected:
    
    bool initialize();
    
    const EMatrix& F1Matrix(const State& state, double  length);
    
    const EMatrix& dr_dx0(const State& state, double length);
    
    const EMatrix& dx_ds(const State& state, double length);
    
  private:
    
    bool parameters();
    
    HelixEquation& helix()
      {return tool<HelixEquation>("equation");}
    
  private:
    
    EMatrix _F1;

    EMatrix _dr_dx0;

    EMatrix _dx_ds;

    int _sdim;

  };
  
}
#endif


