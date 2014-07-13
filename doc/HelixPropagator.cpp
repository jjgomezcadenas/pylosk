#include <recpack/HelixPropagator.h>

//*******************************************
const EMatrix& Recpack::HelixPropagator::F1Matrix(const State& state,
						  double length){
//*******************************************

  _length = length;

  parameters();

	//should be of default representation
  size_t ndim = state.dim();
  _F1 = EMatrix(ndim,ndim,1);

  HelixEquation& param = helix();

  double Bmod       = param.Bmod();

  //--- B=0 case ----
  if (EGeo::is_zero_Bfield(Bmod)){
    EMatrix dr_du0 = EMatrix(3,3,1)*length;
    set_box(dr_du0,0,3,_F1);

    if (verbosity(VVERBOSE)) {
      std::cout << spc(8) << "HelixPropagator::F1Matrix(). For B=0" << std::endl; 
      std::cout << spc(8) << "  F1 \t" << _F1 << std::endl;
      std::cout << spc(8) << "  dr_du0 \t" << dr_du0 << std::endl;
    }

    return _F1;
  }


  double QQ       = param.QQ();
  double theta    = param.theta();
  double sintheta = sin(theta);
  double costheta = cos(theta);

  double qop         = param.q_over_p();
  const EVector& b   = param.b();
  const EVector& u   = param.u();
  const EVector& r0  = param.r0();
  const EVector& r   = param.r();
  const EVector& N   = param.N();

  double A =  (theta-sintheta)/QQ;
  double B =  sintheta/QQ;
  double C = (1-costheta)/QQ;
  double Ap = 1-costheta;
  double Bp = costheta;
  double Cp = sintheta;

  EMatrix dr_du0(3,3,0);     // before was  dr_du0(2,3,0);
  EMatrix du_du0(3,3,0);
  //  EMatrix du_dqop0(3,1,0);

  //---------------------------------
  //  EMatrix bb(b);
  //  EMatrix dr_du0_bis = A*bb*transpose(bb) + B*EMatrix(6,6,1) + C*;

  dr_du0[0][0] = A*b[0]*b[0] + B;
  dr_du0[0][1] = A*b[0]*b[1] - C*b[2];
  dr_du0[0][2] = A*b[0]*b[2] + C*b[1];

  dr_du0[1][0] = A*b[1]*b[0] + C*b[2];
  dr_du0[1][1] = A*b[1]*b[1] + B;
  dr_du0[1][2] = A*b[1]*b[2] - C*b[0];
    
  dr_du0[2][0] = A*b[2]*b[0] - C*b[1];
  dr_du0[2][1] = A*b[2]*b[1] + C*b[0];
  dr_du0[2][2] = A*b[2]*b[2] + B;

  //  std::cout << dr_du0_bis << std::endl;
  //  std::cout << dr_du0 << std::endl;

  //---------------------------------
  du_du0[0][0] = Ap*b[0]*b[0] + Bp;
  du_du0[0][1] = Ap*b[0]*b[1] - Cp*b[2];
  du_du0[0][2] = Ap*b[0]*b[2] + Cp*b[1];
    
  du_du0[1][0] = Ap*b[1]*b[0] + Cp*b[2];
  du_du0[1][1] = Ap*b[1]*b[1] + Bp;
  du_du0[1][2] = Ap*b[1]*b[2] - Cp*b[0];

  du_du0[2][0] = Ap*b[2]*b[0] - Cp*b[1];
  du_du0[2][1] = Ap*b[2]*b[1] + Cp*b[0];
  du_du0[2][2] = Ap*b[2]*b[2] + Bp;
  
  //----------------------------------

  EVector du_dqop0v = theta*N/qop;
  EMatrix du_dqop0(du_dqop0v);

  //----------------------------------


  EVector dr_dqop0v = (r0-r + length*u)/qop;
  EMatrix dr_dqop0(dr_dqop0v);

  set_box(dr_dqop0,  0, 6, _F1);
  set_box(du_dqop0,  3, 6, _F1);

  set_box(dr_du0,   0, 3, _F1);
  set_box(du_du0,   3, 3, _F1);
  
  if (verbosity(VVERBOSE)) {
    std::cout << spc(8) << "HelixPropagator::F1Matrix()" << std::endl; 
    std::cout << spc(8) << "  F1 \t" << _F1 << std::endl;
    std::cout << spc(8) << "  dr_du0 \t" << dr_du0 << std::endl;
    std::cout << spc(8) << "  du_du0 \t" << du_du0 << std::endl;
    std::cout << spc(8) << "  dr_dqop0 \t" << dr_dqop0 << std::endl;
    std::cout << spc(8) << "  du_dqop0 \t" << du_dqop0 << std::endl;
  }

  return _F1;    
}

//*******************************************
const EMatrix& Recpack::HelixPropagator::dx_ds(const State& state, 
					       double length) {
//*******************************************

  _length = length;
  parameters();
  size_t ndim = state.dim();
  _dx_ds = EMatrix(ndim,1,0);
  HelixEquation& param = helix();

  double Bmod = param.Bmod();

  //--- B=0 case ----
  if (EGeo::is_zero_Bfield(Bmod)){
    const EVector& u = param.u();
    for (size_t i = 0; i < 3; i++) _dx_ds[i][0] = u[i];
    
    if (verbosity(VVERBOSE)) 
      std::cout << "HelixPropagator::dx_ds(). \t" << _dx_ds << std::endl;

    return _dx_ds;    
  }

  const EVector& u = param.u();
  const EVector& N = param.N();
  double QQ = param.QQ();

  for (size_t i = 0; i < 3; i++){ 
    _dx_ds[i][0] = u[i];
    _dx_ds[i+3][0] = QQ*N[i];
  }

  if (verbosity(VVERBOSE)) 
      std::cout << spc(8) << "HelixPropagator::dx_ds(). \t" << _dx_ds << std::endl;
  
  return _dx_ds;
}

//*******************************************
const EMatrix& Recpack::HelixPropagator::dr_dx0(const State& state, 
						double length) {
//*******************************************

	size_t ndim = state.dim();
  F1Matrix(state,length);
  _dr_dx0 = box(_F1,0,0,3,ndim);
  if (verbosity(VVERBOSE)) 
    std::cout << spc(8) << "HelixPropagator::dr_dx0() \t" << _dr_dx0 << std::endl;
  return _dr_dx0;
}

//*******************************************
bool Recpack::HelixPropagator::parameters() {
//*******************************************

  return helix().compute_parameters(*_state,_length);
}

