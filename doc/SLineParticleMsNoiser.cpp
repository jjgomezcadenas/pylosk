#include <recpack/SLineParticleMsNoiser.h>
#include <recpack/RayTool.h>
#include <recpack/RepManager.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <recpack/Definitions.h>

//*******************************************
Recpack::SLineParticleMsNoiser::SLineParticleMsNoiser(double x0, size_t ipos) : MsNoiser(x0){
	//*******************************************

	_i = ipos;
	use_tool("equation");
	use_tool("representation");
}

//*******************************************
const EMatrix& Recpack::SLineParticleMsNoiser::QMatrix(const State& state, 
		double length) {
	//*******************************************

	//! compute and return the Q matrix

	// get the running component
	EVector u = equation().direction(state,length);
	unsigned int run = RayTool::running_component(u);

	// get the representation
	Key rep = RayTool::representation(RP::slopes,run);

	//size_t dim = 7;
	//  size_t rdim = representation(run).dim();
	size_t dim = Recpack::rep_default().rep_dim();

	size_t rdim = 5;
	_Q = EMatrix(dim,dim,0);
	EMatrix Q(rdim,rdim,0);

	int index_qoverp;
	bool ok = Recpack::rep_default().index(RP::qoverp, index_qoverp);
	if (!ok){
		if (verbosity(VVERBOSE)){
			std::cout << " SLineParticleMsNoiser::QMatrix: no qoverp in default representation: "<<
				Recpack::rep().default_rep_name()<< std::endl;
		}
		return _Q;
	}

	// retrieve the total momentum
	double qoverp = fabs(state.vector()[index_qoverp]);

	if (qoverp == 0. || fabs(_x0 == 0.)){
		if (verbosity(NORMAL)){
			std::cout << " SLineParticleMsNoiser::QMatrix: no momentum or no x0 in noiser" 
				<<" q/p, X0 = " << qoverp <<  _x0 << std::endl;
		}
		return _Q;
	}



	size_t sdim = u.num_row();
	EVector rp = RayTool::rp(run,u);

	// retrieve lambda
	double lambda0 = RayTool::lambda(run,u);

	double p = fabs(1/qoverp);
	double x0 = _x0;
  
  Key PID= "";
	if (state.names().has_key(RP::PID)) 
	  PID = state.name(RP::PID);

  double m=0;
  if (PID=="Muon")
	  m = _mmuon;
	else if (PID=="Electron")
	  m = _melectron;
	else if (PID=="Pion")
		m = _mpion;
	else if (PID=="Proton")
		m = _mproton;
	else if (PID=="Kaon")
		m = _mkaon;
	else
		m = _mmuon;

	double beta2 = p*p/(p*p+m*m);
  double beta = 1.;
  if(beta2>0.){
    beta = sqrt(beta2);
  }
    
	
	double k=pow((13.6*MeV/(p*beta)),2)/x0;
  

	double ll = fabs(length);

	for (size_t j = 0; j < sdim-1; j++) {

		unsigned int ind = RayTool::rp_index(run,j);
		double cos4_theta0=pow(1./(1+rp[j]*rp[j]),2);
		double lambda0_2 = lambda0*lambda0;

		// position components are always x,y,z. Slope component may change
		Q[ind][ind] =                     k*lambda0_2*(pow(ll,3))/(3*cos4_theta0);
		Q[j+_i][ind] = Q[ind][j+_i] = k*lambda0*(pow(ll,2))/(2*cos4_theta0);
		Q[j+_i][j+_i] =               k*ll/cos4_theta0;

	}

	// gets the jacobian slopes:z --> default
	//EMatrix J = representation().jacobian_inverse(state,rep);  

	// gets the jacobian slopes_curv --> current default

	EMatrix J = Recpack::rep().jacobian_inverse(state,rep);     



	_Q = J*Q*transpose(J);

	if (verbosity(VVERBOSE)) {
		std::cout << "  --------------------------------" << std::endl;
		std::cout << "  SLineParticleMsNoiser::QMatrix  computing ms : "
			<< " X0 = " << x0 << ", length = " << length 
			<< ", p = " << p << ", PID = " << PID << ",  beta"<< beta <<", run = " << run << std::endl;
		std::cout << "  Q matrix in slopes representation: " << Q << std::endl;
		std::cout << "  Jacobian to default representation: " << J << std::endl;
		std::cout << "  Q matrix: " << _Q << std::endl;
		std::cout << "  --------------------------------" << std::endl;
	}
	else if (verbosity(VERBOSE)) {
		std::cout << "  --------------------------------" << std::endl;
		std::cout << "  SLineParticleMsNoiser::QMatrix  computing ms : "
			<< " X0 = " << x0 << ", length = " << length 
			<< ", p = " << p << ", PID " << PID << ", beta " << beta << std::endl;
    std::cout << "  Q matrix: " << _Q << std::endl;
		std::cout << "  --------------------------------" << std::endl;
	}
	else if (verbosity(DETAILED)) {
		std::cout << "    SLineParticleMsNoiser::QMatrix(). "
			<< "X0 = " << x0 << ", l = " << length
			<< ", p = " << p << ", PID " << PID << ", beta " << beta
			<< ", Q = " << _Q[0][0] << std::endl;
	}  

	return _Q;
}

