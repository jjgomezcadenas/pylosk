#ifndef SLinePropagator_h
#define SLinePropagator_h 1

#include <recpack/LengthPropagator.h>
#include <recpack/SLineEquation.h>

using namespace Recpack;

namespace Recpack {

  //! propagator for the stragiht line
  class SLinePropagator : public LengthPropagator {
  public:

    //! default constructor
    SLinePropagator(size_t sdim) 
      {_sdim = sdim;}
    
    //! default destructor
    virtual ~SLinePropagator() {}
    
  protected:
  
    const EMatrix& F1Matrix(const State& state, double  length);
    
    const EMatrix& dr_dx0(const State& state, double length);

    const EMatrix& dx_ds(const State& state, double length);

  private:

    SLineEquation& equation()
      {return tool<SLineEquation>("equation");}

  protected:

    size_t _sdim;

    EMatrix _F1;

    EMatrix _dr_dx0;

    EMatrix _dx_ds;

  };

}

#endif


