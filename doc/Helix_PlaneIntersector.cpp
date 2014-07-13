#include <recpack/Helix_PlaneIntersector.h>
#include <CLHEP/Units/PhysicalConstants.h>

//*******************************************
Recpack::Helix_PlaneIntersector::Helix_PlaneIntersector() {
//*******************************************

  use_tool("equation"); _helix = NULL;
  use_tool("intersector/numerical");
  use_tool("intersector/helix_num");
  use_tool("intersector/sline/plane");
  
}

//*******************************************
bool Recpack::Helix_PlaneIntersector::length(const State& state,
					     const Surface& surf,
					     double& length){
//*******************************************

  bool ok = initialize(state,surf);
  if (!ok){
    if (verbosity(DETAILED))
      std::cout << " Helix_PlaneIntersector::length(): Failed in initialization !!!" << std::endl;
    return false;
  }
  length = 0;

  double Bmod = _helix->Bmod();

  //------- B=0 case --------
  if (EGeo::is_zero_Bfield(Bmod)){
    if (verbosity(DETAILED)){
      std::cout << " Helix_PlaneIntersector::length(): B=0 --> Use SLine_PlaneIntersector" << std::endl;
    }
    ok = sline().length(state,surf,length);

    // we need to copy ds_dS from the SLine intersector to the Helix intersector
    _ds_dS = sline().ds_dS(state,surf,length);
  }
  //------- normal B case --------
  else{
    double proj = fabs(dot(_np,_b));
    
    if (verbosity(VERBOSE)){ 
      std::cout << spc(8) << "Helix_PlaneIntersector::length(). dot(b,n) = " << proj << std::endl;
    }
    
    if ( proj > 1-EGeo::zero_transverse())  ok = caseBParallel(length);
    else if ( proj < EGeo::zero_transverse()) ok = caseBPerpendicular(length);
    else{
      if (verbosity(VVERBOSE)){ 
	std::cout << spc(8) << "  --> Plane normal and B field are not parallel nor perpendicular. Intersect numerically " << std::endl;
      }
      ok = numerical().length(state,surf,length);
    }
  }


  if (verbosity(VERBOSE)){ 
    if (!ok) std::cout << spc(8) << "Helix_PlaneIntersector::length()  "
		       << " can not compute path !!!!"  << std::endl;
    else     std::cout << spc(8) << "Helix_PlaneIntersector::length(). "  
		       << "Intersection OK. length = " << length << std::endl;
  }

  if (!ok) return false;

  const EVector& pos = _helix->position(state,length);
  const EVector& dir = _helix->direction(state,length);

  if (verbosity(VVERBOSE))
    std::cout << spc(8) << "Helix_PlaneIntersector::length(): extrapolated pos = " << print(pos) << ", dir = " << print(dir) << std::endl;

  bool inside = surf.is_inside(pos, verbosity());

  return inside;
}

//*******************************************
bool Recpack::Helix_PlaneIntersector::caseBParallel(double& length) {
//*******************************************

  // when the field is parallel to the plane normal (so it is normal to the plane), 
  // there is only one solution for the intersection since in this projection the helix is a sinusoidal 
  // (it will intersect the plane in only one point)

  length = 0;
  if (_gamma == 0.) return false;
  length = dot(_b,_rp-_r0)/_gamma;

  int sdim = _b.num_row();
  EMatrix BB(1,sdim);
  for (int i = 0; i < sdim; i++) BB[0][i] = _b[i];
  
  _ds_dS = BB/_gamma;

  if (verbosity(VVERBOSE)) {
    std::cout << spc(8) << "Helix_PlaneIntersector::caseBParallel(). Allow zero length = " << _allow_zero_length << ", l sign = " << _length_sign << std::endl;
    std::cout << spc(8) << "  length = " << length   << std::endl;
    std::cout << spc(8) << "  b      = " << print(_b) << std::endl;
    std::cout << spc(8) << "  gama   = " << _gamma   << std::endl;
    std::cout << spc(8) << "  rp     = " << print(_rp) << std::endl;
    std::cout << spc(8) << "  r0     = " << print(_r0) << std::endl;
    std::cout << spc(8) << "  ds_dS  = " <<_ds_dS << std::endl;
  }

  if (EGeo::is_zero_l(length) && !_allow_zero_length) return false;
  else if (length>0 && _length_sign==-1) return false;
  else if (length<0 && _length_sign==1) return false;

  return true;
}

//*******************************************
bool Recpack::Helix_PlaneIntersector::caseBPerpendicular(double& length) {
//*******************************************

  // when the field is normal to the plane normal (so it is parallel to the plane), 
  // there are two solutions for the intersection since in this projection the helix is a circle 
  // (it will intersect the plane in two points)

  bool ok=false;
  length = 0;

  // angle between n0 (normal to trajectory) and normal to plane
  long double cos_theta0 = dot(_n0,_np);
  long double theta0 = acos(cos_theta0);

  // resolve the ambigutity in theta
  long double sin_theta0 = -dot(_m0,_np);
  if (sin_theta0 < 0) theta0 *=-1;

  // This is the intersection between the helix and the plane
  long double cos_theta = -dot(_np,_rp-_r0)/_RR + cos_theta0;
  long double theta = 0;

  // compute sin(theta)
  EVector np_perp = crossprod(_np,_b);
  long double sin_theta = -dot(np_perp,_rp-_r0)/_RR + sin_theta0;

  long double vtheta[4];
  long double vlength[4];
  long double lpos_min=1e10;
  long double lneg_min=1e10;
  int ipos_min=-1;
  int ineg_min=-1;
  long double lpos_max=0.;
  long double lneg_max=0.;
  int ipos_max=-1;
  int ineg_max=-1;

  // If |cos_theta|>1 there is no solution --> the surface is not intersected
  if (fabs(cos_theta)<1){
    theta = acos(cos_theta); 

    // The  first solution (1st point in space)
    vtheta[0] = theta;

    // The second solution (2nd point in space) would be symmetric (with -sin_theta)
    vtheta[1] = 2*CLHEP::pi-theta;

    // Same point in space as 2nd solution but differen theta
    vtheta[2] = -vtheta[0];

    // Same point in space as 1st solution but differen theta
    vtheta[3] = -vtheta[1];

    // Find the min and max length solutions for positive and negative length
    for (int i=0;i<4;i++){
      vlength[i] = -1/_QQ*(vtheta[i] - theta0);
      if (EGeo::is_zero_l(vlength[i]) && !_allow_zero_length) continue;
      if (fabs(vtheta[i]-theta0) > 2*CLHEP::pi) continue;  // solutions above 2pi are discarded 
      if (vlength[i]>=0 || EGeo::is_zero_l(vlength[i])){
	if (fabs(vlength[i])<lpos_min){
	  lpos_min= fabs(vlength[i]);
	  ipos_min = i;
	}
	if (fabs(vlength[i])>=lpos_max){
	  lpos_max= fabs(vlength[i]);
	  ipos_max = i;
	}
      }
      if (vlength[i]<=0 || EGeo::is_zero_l(vlength[i])){
	if (fabs(vlength[i])<lneg_min){
	  lneg_min= fabs(vlength[i]);
	  ineg_min = i;
	}
	if (fabs(vlength[i])>=lneg_max){
	  lneg_max= fabs(vlength[i]);
	  ineg_max = i;
	}
      }
    }

    // For l sign 2 and 3 allow a minimum difference in angle (this is for curv back tracks). 
    // Harcoded for the moment
    double dtheta_min = 0;
    if (abs(_length_sign)>1) dtheta_min = CLHEP::pi/2.2;

    // Choose the right solution
    long double lpos;
    long double lneg;
    
    if (fabs(vtheta[ipos_min]-theta0) >= dtheta_min)
      lpos = vlength[ipos_min];
    else
      lpos = vlength[ipos_max];
    if (fabs(vtheta[ineg_min]-theta0) >= dtheta_min)
      lneg = vlength[ineg_min];
    else
      lneg = vlength[ineg_max];


    if (_length_sign==1 || _length_sign==2)
      length = lpos;
    else if (_length_sign==-1 || _length_sign==-2)
      length = lneg;
    else{
      if (fabs(lneg) < fabs(lpos))
	length = lneg;
      else 
	length = lpos;
    }
    ok = true;

    if (ok){
      int sdim = _np.num_row();
      EMatrix NP(1,sdim);
      for (int i = 0; i < sdim; i++) NP[0][i] = _np[i];
      
      _ds_dS = 1/(_QQ*sqrt(1-cos_theta*cos_theta)*_RR)*NP;
    }

  }

  if (verbosity(VVERBOSE)) {
    std::cout << spc(8) << "Helix_PlaneIntersector::caseBPerpendicular(). Allow zero length = " << _allow_zero_length << ", l sign = " << _length_sign << std::endl;
    std::cout << spc(8) << "  length          = " << length << std::endl;
    std::cout << spc(8) << "  length1         = " << vlength[0] << std::endl;
    std::cout << spc(8) << "  length2         = " << vlength[1] << std::endl;
    std::cout << spc(8) << "  length3         = " << vlength[2] << std::endl;
    std::cout << spc(8) << "  length4         = " << vlength[3] << std::endl;
    std::cout << spc(8) << "  cos(theta0)     = " << cos_theta0 << std::endl;
    std::cout << spc(8) << "  cos(theta)      = " << cos_theta << std::endl;
    std::cout << spc(8) << "  sin(theta0)     = " << sin_theta0 << std::endl;
    std::cout << spc(8) << "  sin(theta)      = " << sin_theta << std::endl;
    std::cout << spc(8) << "  sin(theta0)'    = " << sqrt(1-cos_theta0*cos_theta0) << std::endl;
    std::cout << spc(8) << "  sin(theta)'     = " << sqrt(1-cos_theta*cos_theta) << std::endl;
    std::cout << spc(8) << "  theta0          = " << theta0 << std::endl;
    std::cout << spc(8) << "  theta           = " << theta << std::endl;
    std::cout << spc(8) << "  theta1          = " << vtheta[0] << std::endl;
    std::cout << spc(8) << "  theta2          = " << vtheta[1] << std::endl;
    std::cout << spc(8) << "  theta3          = " << vtheta[2] << std::endl;
    std::cout << spc(8) << "  theta4          = " << vtheta[3] << std::endl;
    std::cout << spc(8) << "  theta-theta0    = " << theta-theta0 << std::endl;
    std::cout << spc(8) << "  dot(np,rp-r0)/R = " << dot(_np,_rp-_r0)/_RR << std::endl;
    std::cout << spc(8) << "  dot(np_p,rp-r0)/R = " << dot(np_perp,_rp-_r0)/_RR << std::endl;
    std::cout << spc(8) << "  np              = " << print(_np) << std::endl;
    std::cout << spc(8) << "  np_perp         = " << print(np_perp) << std::endl;
    std::cout << spc(8) << "  rp              = " << print(_rp) << std::endl;
    std::cout << spc(8) << "  r0              = " << print(_r0) << std::endl;
    std::cout << spc(8) << "  u0              = " << print(_u0) << std::endl;
    std::cout << spc(8) << "  n0              = " << print(_n0) << std::endl;
    std::cout << spc(8) << "  m0              = " << print(_m0) << std::endl;
    std::cout << spc(8) << "  R               = " << _RR << std::endl;
    std::cout << spc(8) << "  Q               = " << _QQ << std::endl;
    std::cout << spc(8) << "  charge          = " << _charge << std::endl;
    std::cout << spc(8) << "  ds_dS           = " << _ds_dS << std::endl;
  }

  if (!ok){
    if (verbosity(VERBOSE)) 
      std::cout << spc(8) << "--WARNING: HelixPlane::CaseBPerperdicular No intersection" 
		<< std::endl;
    return false;      
  } 

  return true;
}

//*******************************************
const EMatrix& Recpack::Helix_PlaneIntersector::ds_dS(const State& state, const Surface& surf, double length){
//*******************************************

  if (verbosity(VVERBOSE)) {
    std::cout << spc(8) << "Helix_PlaneIntersector::ds_dS()" << std::endl;
    std::cout << spc(8) << "  ds_dS          = " << _ds_dS << std::endl;
  }

  return _ds_dS;

}

//*******************************************
bool Recpack::Helix_PlaneIntersector::initialize(const State& state, 
						 const Surface& surf) {
//*******************************************

  if (_helix == NULL) 
    _helix = &(tool<HelixEquation>("equation"));
 
  _rp = surf.position();
  _np = surf.normal();

  bool ok =_helix->compute_parameters(state,0.);
  if (!ok) return false;

  _b  = _helix->b();
  _r0 = _helix->r0();
  _u0 = _helix->u0();
  _n0 = _helix->n0();
  _m0 = _helix->m0();

  _alpha = _helix->alpha();
  _gamma = _helix->gamma();
  _QQ = _helix->QQ();
  _RR = _helix->RR();
  _charge = _helix->charge();

  return true;
}
