#ifndef LengthPropagator_h
#define LengthPropagator_h 1

#include <recpack/Propagator.h>

using namespace Recpack;

namespace Recpack {

  //! propagator a state a given path length
  class LengthPropagator : public Propagator {    
  public:
    
    //! default constructor
    LengthPropagator() {}
    
    //! default destructor
    virtual ~LengthPropagator() {} 

    //! Matrix (partial derivative dr/dx0 r=geometrical position x0=at origin)
    virtual const EMatrix& dr_dx0(const State& state, double length) = 0;

  protected:

    //! Matrix (partial derivative dx/dx0 x=vector parameters x0= at origin)
    //! it depends on the propagation model
    virtual const EMatrix& F1Matrix(const State& state, double length) = 0;

    
    //! Matrix (partial derivative dx/ds x=vector parameters, s=path_length)
    virtual const EMatrix& dx_ds(const State& state, double length) = 0;

    //! F Matrix (absolute derivative dx/dx0 x=vector parameters x0 at origin)
    //! length parameters delegate F = F1 + Fs = (dx/dx0+dx/dx0*dx0/ds)
    //! where s is the path length. Fs is computed if propagtion to a surface
    const EMatrix& transport_matrix(const State& state, double length);
    
    
    //! Matrix (dx/dx0*dx0/dx x=vector parameters x0= at origin)
    //! note: only computed if the propagation is to a surface
    const EMatrix& FsMatrix(const State& state, const Surface& surface,
			    double length);


    //! Matrix (partial derivative dx/dS x=vector parameters, s=surface parameters)
    const EMatrix& dx_dS(const State& state, const Surface& surface, 
    			 double length);    
  protected:
    
    //! Fs matrix
    EMatrix _Fs;

  };

}
#endif




