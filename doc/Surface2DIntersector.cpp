#include <recpack/Surface2DIntersector.h>


//*******************************************
Recpack::Surface2DIntersector::Surface2DIntersector(){
//*******************************************

  use_tool("equation"); 
  use_tool("propagator"); 

}

//*******************************************
const EMatrix& Recpack::Surface2DIntersector::ds_dx0(const State& state, 
						     const Surface& surface,
						     double length){
//*******************************************


  const EVector& u = equation().direction(state,length);
  const EVector& r = equation().position(state,length);
  const EVector& np = surface.normal(r);
  int sdim = np.num_row();
  EMatrix NP(1,sdim);
  for (int i = 0; i < sdim; i++) NP[0][i] = np[i];
  double deno = dot(np,u);
  if (deno == 0){
    if (verbosity(DETAILED))  
      std::cout << "    Surface2DIntersector::ds_dx0(). Can not compute ds_dx0. Surface and direction vector are parallel." << std::endl;
    size_t ndim = state.dim();
    _ds_dx0 = EMatrix(1,ndim,0);
  }
  else{
    const EMatrix& drdx0 = propagator().dr_dx0(state,length);
    // testing it was [-] (NP*)
    _ds_dx0 =  - (NP*drdx0)/(deno);
  }

  if (verbosity(VVERBOSE)) 
    std::cout << "    Surface2DIntersector::ds_dx0(). ds_dx0  = " << _ds_dx0 << std::endl;
  return _ds_dx0;

}
