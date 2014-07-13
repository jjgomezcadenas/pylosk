#include <recpack/LengthPropagator.h>


//*******************************************
const EMatrix& Recpack::LengthPropagator::transport_matrix(const State& state, 
							   double length) {
//*******************************************

  _F = F1Matrix(state,length);
  if (_status["via_surface"]) 
    _F += FsMatrix(state,*_surface, length);
  if (verbosity(VERBOSE)){
    std::cout << "    LengthPropagator::transport_matrix(). F matrix  = " << std::endl;
    std::cout << _F << std::endl;
  }

  // TODO
  for (int i=0;i<_F.num_row();i++)
    for (int j=0;j<_F.num_row();j++)
      if (fabs(_F[i][j]) < 1e-10) _F[i][j] = 0;
      
  return _F;
}

//*******************************************
const EMatrix& Recpack::LengthPropagator::FsMatrix(const State& state, 
						   const Surface& surface, 
						   double length) {
//*******************************************

  _Fs = dx_ds(state,length) * intersector().ds_dx0(state,surface,length);
  if (verbosity(VVERBOSE)){ 
    std::cout << "    LengthPropagator::FsMatrix(). Fs  = " << std::endl;
    std::cout << _Fs << std::endl;
  }
  return _Fs;
}


//*******************************************
const EMatrix& Recpack::LengthPropagator::dx_dS(const State& state, 
						const Surface& surface, 
						double length){
							 
//*******************************************

  _dx_dS = dx_ds(state,length) * intersector().ds_dS(state,surface,length);
  if (verbosity(VVERBOSE)){ 
    std::cout << "    LengthPropagator::dx_dS(). dx_dS  = " << std::endl;
    std::cout << _dx_dS << std::endl;
  }
  return _dx_dS;
}
