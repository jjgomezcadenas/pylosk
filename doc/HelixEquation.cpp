#include <recpack/HelixEquation.h>
#include <recpack/RepManager.h>
#include <recpack/EGeo.h>
#include <recpack/Definitions.h>

//*******************************************
Recpack::HelixEquation::HelixEquation() {
	//*******************************************

	size_t dim = Recpack::rep_default().rep_dim();

	_B = EVector(3,0); _x0 = EVector(dim,0); _length = EGeo::inf_l();

	_status[RP::ready]=false;

	use_tool("representation"); 

}

//*******************************************
bool Recpack::HelixEquation::compute_parameters(const State& cstate,
		double length) {
	//*******************************************

	// make a copy of the initial state. Otherwise it cannot be modified by the representation
	State state = cstate;

	//if (state.name(RP::representation) != RP::default_rep)
  //representation().convert(state,RP::default_rep);
	
	// converts to the default representation
	if (state.name(RP::representation) != Recpack::rep().default_rep_name())
		RP::rep().convert(state,Recpack::rep().default_rep_name());
	

	const EVector& vec = state.vector();
	if (vec == _x && EGeo::is_zero_l(length-_length)){
		return true;
	}  

	if (!_status[RP::ready]){
		if (verbosity(DETAILED)) {
			std::cout << spc(8) << "HelixEquation::compute_parameters(). No property BFiels is available in this volume" << std::endl;     
		}
		return false;
	}

	_x0 = vec;
	_length = length;

	// read the parameters from the vector
	bool ok = read_vector(state);
	if (!ok) return false;


	_Bmod = _B.norm();
	/*
		 if (_Bmod == 0.) {
		 msg(WARNING," HelixEquation:  null magnetic field! ");
		 return false;
		 }
		 */
 
	//compute time <-- needs p and length
  Key PID= "";
  if (state.names().has_key(RP::PID)) 
  	PID = state.name(RP::PID);


  double mass=0;
  if (PID=="Muon")
  	mass = _mmuon;
  else if (PID=="Electron")
    mass = _melectron;
  else if (PID=="Pion")
    mass = _mpion;
  else if (PID=="Proton")
    mass = _mproton;
  else if (PID=="Kaon")
    mass = _mkaon;
  else
    mass = _mmuon;

  
 	double beta = _p0/sqrt(_p0*_p0+mass*mass);

  double v_particle = beta*_speed_light; 
    
  if (v_particle!=0)
    		_t = _t0 + _length/v_particle;
  
	//----  B=0 case ----//
	if (EGeo::is_zero_Bfield(_Bmod)){
		// compute position
		_r = _r0 + _length*_u0;

		// compute direction
		_u = _u0;        

		// compute vector
		load_vector(state);

		if (verbosity(VERBOSE)) {
			std::cout << spc(8) << "HelixEquation::compute_parameters()" << std::endl; 
			std::cout << spc(8) << "  length = " << _length << std::endl;
			std::cout << spc(8) << "  r      = " << print(_r) <<std::endl;      
		}

		return true;
	}

	_b = _B/_Bmod;

	// natural units: MeV, ktesla, & mm
	// 0.3 [GeV]/([T][m]) = 300
	double factor = 0.3*GeV/(CLHEP::tesla*CLHEP::m);
	_QQ= -_charge*factor*(_Bmod/_p0);

	// cos of the angle between momentum and magnetic field
	_gamma = dot(_b,_u0);
	// sin of the angle between momentum and magnetic field
	_alpha = sqrt(1-_gamma*_gamma);

	// radius of the circle
	if (_QQ !=0) _RR= _alpha/_QQ;
	else _RR=100000;


	_N0 = crossprod(_b,_u0);  

	_n0 = crossprod(_b,_u0);
	if (_alpha != 0) _n0 /= _alpha;


	_m0 = crossprod(_n0,_b);
	double mm0 = _m0.norm();
	if (mm0 != 0) _m0 /= mm0;  

	// compute position
	_theta    = _QQ*length;
	_costheta = cos(_theta);
	_sintheta = sin(_theta); 

	_r = _r0 + _gamma/_QQ*(_theta-_sintheta)*_b 
		+ (1/_QQ)*_sintheta*_u0 + _alpha/_QQ*(1-_costheta)*_n0;

	// compute direction
	_u = _gamma*(1-_costheta)*_b
		+ _costheta*_u0 + _alpha*_sintheta*_n0;
	_N = crossprod(_b,_u);

	// compute vector
	load_vector(state);

	if (verbosity(VVERBOSE)) {
		std::cout << spc(8) << "HelixEquation::compute_parameters()" << std::endl; 
		std::cout << spc(8) << "  B    \t" << print(_B) << std::endl;
		std::cout << spc(8) << "  Bmod \t" << _Bmod << std::endl;
		std::cout << spc(8) << "  r0   \t"   << print(_r0) << std::endl; 
		std::cout << spc(8) << "  u0   \t"   << print(_u0) << std::endl; 
		std::cout << spc(8) << "  n0   \t"    << print(_n0) << std::endl; 
		std::cout << spc(8) << "  m0   \t"    << print(_m0) << std::endl; 
		std::cout << spc(8) << "  qop0  \t"  << _qop0 << std::endl; 
	  std::cout << spc(8) << "  t0  \t"  << _t0 << std::endl;
		std::cout << spc(8) << "  p    \t"   << _p0 << std::endl; 
		std::cout << spc(8) << "  charge\t"<< _charge << std::endl; 
		std::cout << spc(8) << "  QQ   \t"    << _QQ << std::endl; 
		std::cout << spc(8) << "  RR   \t"    << _RR << std::endl; 
		std::cout << spc(8) << "  gamma\t" << _gamma << std::endl; 
		std::cout << spc(8) << "  alpha\t"  << _alpha << std::endl; 
		std::cout << spc(8) << "  at length : " << _length << std::endl;
		std::cout << spc(8) << "  theta\t"  << _theta << std::endl; 
		std::cout << spc(8) << "  costheta\t"  << _costheta << std::endl; 
		std::cout << spc(8) << "  sintheta\t"  << _sintheta << std::endl; 
		std::cout << spc(8) << "  r    \t"  << print(_r) << std::endl;
		std::cout << spc(8) << "  t  \t"  << _t << std::endl; 
		std::cout << spc(8) << "  u    \t"  << print(_u) << std::endl; 
		std::cout << spc(8) << "  N    \t"  << print(_N) << std::endl;
		std::cout << spc(8) << "  state : >> " << std::endl;
		std::cout << spc(8) << "  vector \t" << print(_x) << std::endl;
		std::cout << spc(8) << "  par    \t" << _r[2] << std::endl;
	}

	return true;
}


//*******************************************
bool Recpack::HelixEquation::read_vector(const State& cstate) {
	//*******************************************

	// make a copy of the initial state. Otherwise it cannot be modified by the representation
	State state = cstate;

 	//if (state.name(RP::representation) != RP::default_rep)
  //representation().convert(state,RP::default_rep);
	
	// converts to the default representation
	if (state.name(RP::representation) != Recpack::rep().default_rep_name())
		Recpack::rep().convert(state, Recpack::rep().default_rep_name()); 

	_x0 = state.vector();

	// position
	_r0= box(_x0,0,3);

	// direction
	//  _rp0 = box(vec,3,2);
	//  _u0 = RayTool::u(_rp0);
	_u0 = box(_x0,3,3);

	// q_over_p
	_qop0 = box(_x0,6,1)[0];

	//time (if present)
	int index_t;
	bool ok = Recpack::rep_default().index(RP::time, index_t);
	if (ok)
		_t0 = _x0[index_t];
	else
		_t0 = 0.;	

	// take the charge and the 
	if (fabs(_qop0) <= 0.) {
		msg(WARNING,"Helix Equation::read_vector(). null 1/p!");
		return false;
	}

	// assume charge +-1  
	_charge = _qop0/fabs(_qop0);
	_p0 = fabs(1/_qop0);


	if (verbosity(VVERBOSE)) {
		std::cout << spc(8) << "HelixEquation::read_vector()" << std::endl; 
		std::cout << spc(8) << "  x0  = " << print(_x0) << std::endl; 
		std::cout << spc(8) << "  r0  = " << print(_r0) << std::endl; 
		std::cout << spc(8) << "  u0  = " << print(_u0) << std::endl; 
		std::cout << spc(8) << "  qop0 = " << _qop0 << std::endl; 
		std::cout << spc(8) << "  p0  = " << _p0 << std::endl; 
		std::cout << spc(8) << "  t0  \t"  << _t0 << std::endl; 
		std::cout << spc(8) << "  q   = " << _charge << std::endl; 
	}

	return true;

}

//*******************************************
void Recpack::HelixEquation::load_vector(const State& state) {
	//*******************************************

	_x = _x0;
	set_box(_r,0,_x);
	//  set_box(_rp,3,_x);
	set_box(_u,3,_x);
	//  _x[5]=_charge/_p0;
	_x[6]=_charge/_p0;


	//time (if present)
	int index_t;
	bool ok = Recpack::rep_default().index(RP::time, index_t);
	if(ok)
		_x[index_t] = _t;

	if (verbosity(VVERBOSE)) {
		std::cout << spc(8) << "HelixEquation::load_vector(). x = " << print(_x) << std::endl; 
	}

}
