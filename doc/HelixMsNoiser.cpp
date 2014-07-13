#include <recpack/HelixMsNoiser.h>
#include <recpack/RepManager.h>
#include <recpack/SLineParticleMsNoiser.h>
#include <recpack/RayTool.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <recpack/Definitions.h>

//*******************************************
Recpack::HelixMsNoiser::HelixMsNoiser(double x0, size_t ipos) : MsNoiser(x0){
	//*******************************************

	_i = ipos;
	use_tool("equation");
	use_tool("representation");
	use_tool("propagator");
	use_tool("noiser/sline/ms");
}

//*******************************************
const EMatrix& Recpack::HelixMsNoiser::QMatrix(const State& state, 
		double length) {
	//*******************************************

	//! compute and return the Q matrix 
	//size_t ndim = state.dim();
	//size_t ndim =7;

	size_t ndim = Recpack::rep_default().rep_dim();

	_Q = EMatrix(ndim,ndim,0);

	int index_qoverp;
	bool ok = Recpack::rep_default().index(RP::qoverp, index_qoverp);
	if (!ok){
		if (verbosity(VVERBOSE)){
			std::cout << " HelixMsNoiser::QMatrix: no qoverp in default representation: "<<
				Recpack::rep().default_rep_name()<< std::endl;
		}
		return _Q;
	}

	// retrieve the total momentum
	double qoverp = fabs(state.vector()[index_qoverp]);

	if (qoverp == 0. || fabs(_x0 == 0.)){
		if (verbosity(VVERBOSE)){
			std::cout << " HelixMsNoiser::QMatrix: no momentum or no x0 in noiser" 
				<<" q/p, X0 = " << qoverp <<  _x0 << std::endl;
		}
		return _Q;
	}


	double p = fabs(1/qoverp);
	double x0 = _x0;
	double Bmod = equation().Bmod();

	//------- B=0 case --------
	if (EGeo::is_zero_Bfield(Bmod)){
		if (verbosity(DETAILED)){
			std::cout << " HelixMsNoiser::QMatrix(): B=0 --> Use SLineParticleMsNoiser" << std::endl;
		}
		sline().set_x0(_x0);
		_Q = sline().QMatrix(state,length);
	}
	//------- normal B case --------
	else{
    
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


		EMatrix Qx(ndim,ndim,0);
		EMatrix Qy(ndim,ndim,0);

		//--------------------

		// get the direction vector
		EVector u = equation().direction(state);

		// get the running component
		unsigned int run = RayTool::running_component(u);

		EVector rp = RayTool::rp(run,u);
		double cos4_thetax0=pow(1./(1+rp[0]*rp[0]),2);
		double cos4_thetay0=pow(1./(1+rp[1]*rp[1]),2);


		// gets the jacobian slopes_curv --> default

		Key rep = RayTool::representation(RP::slopes_curv, run);
		//EMatrix J = representation().jacobian_inverse(state,rep);      

		// gets the jacobian slopes_curv --> current default

		EMatrix J = Recpack::rep().jacobian_inverse(state,rep);     

		EMatrix Jx_3(ndim,ndim,0);
		EMatrix Jy_3(ndim,ndim,0);

		for (int i=0;i<7;i++){
			for (int j=0;j<7;j++){
				Jx_3[i][j] = J[i][3]*J[j][3];
				Jy_3[i][j] = J[i][4]*J[j][4];
			}      
		}

		EMatrix Qx_0 = k/cos4_thetax0*fabs(length)*Jx_3; 
		EMatrix Qy_0 = k/cos4_thetay0*fabs(length)*Jy_3; 

		const EMatrix& F = propagator().FMatrix();  
		EMatrix F_T = transpose(F);

		Qx = F*(Qx_0)*F_T;
		Qy = F*(Qy_0)*F_T;


		//---------- Q matrix

		_Q = Qx+Qy;

		if (verbosity(VVERBOSE)) {
			std::cout << "  --------------------------------" << std::endl;
			std::cout << "  HelixMsNoiser::QMatrix  computing ms : "
				<< " X0 = " << x0 << ", length = " << length
				<< ", p = " << p << ", PID = " << PID << ",  beta "<< beta <<", run = " << run << std::endl;
			std::cout << "  Qx matrix: " << Qx << std::endl;
			std::cout << "  Qy matrix: " << Qy << std::endl;
			std::cout << "  Q  matrix: " << _Q << std::endl;
			std::cout << "  --------------------------------" << std::endl;
		}

		else if (verbosity(VERBOSE)) {
			std::cout << "  HelixMsNoiser::QMatrix  computing ms : "
				<< " X0 = " << x0 << ", length = " << length 
				<< ", p = " << p << ", PID " << PID << ", beta " << beta << std::endl;
			std::cout << "  Q matrix: " << _Q << std::endl;
		}
		else if (verbosity(DETAILED)) {
			std::cout << "    HelixMsNoiser::QMatrix(). "
				<< "X0 = " << x0 << ", l = " << length 
				<< ", p = " << p << ", PID " << PID << ", beta " << beta
				<< ", Q = " << _Q[0][0] << std::endl;
		}  
	}

	return _Q;


}

