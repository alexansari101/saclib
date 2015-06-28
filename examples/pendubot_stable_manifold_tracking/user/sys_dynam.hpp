#ifndef SYS_DYNAM_HPP
#define SYS_DYNAM_HPP

namespace sac {

  const size_t xlen = 4;
  const size_t ulen = 1;
  const double p1 = 0.0308, p2 = 0.0106, p3 = 0.0095,
    p4 = 0.2087, p5 = 0.0629;

  std::vector< std::vector< double > > get_csv( const char * filename );
  double s1_intp(double th1, double th2);
  double s2_intp(double th1, double th2);
  mat_type Gxdes(double th1, double th2);


  //[ The rhs of x' = f(x) defined as a class
  // USER SPECIFIED:
  class sys_dynam {
    double g_;
    b_control & u_;
    state_type u_curr_;
  
  public:
    sys_dynam( b_control & u ) : g_(9.81) , u_(u) , u_curr_(ulen) { }
  
    void operator() (const state_type &x, state_type &dxdt, const double t)
    {
      u_(t, u_curr_);
      dxdt[0] = x[1];
      dxdt[1] = ( g_*(2*p2*p4 - p3*p5)*sin(x[0]) + 
		  g_*p3*p5*sin(x[0]-2*x[2]) + 2*p2*u_curr_[0] - 
		  2*p3*sin(x[0]-x[2])*(p3*cos(x[0]-x[2])*pow(x[1], 2) + 
				       p2*pow(x[3], 2))
		  )/( 2*p1*p2 - 2*pow(p3, 2)*pow(cos(x[0]-x[2]), 2) );
      dxdt[2] = x[3];
      dxdt[3] = (-g_*p3*p4*cos(x[0]-x[2])*sin(x[0]) + g_*p1*p5*sin(x[2]) -
		 p3*cos(x[0]-x[2])*u_curr_[0] + p3*sin(x[0]-x[2])*
		 (p1*pow(x[1], 2) + p3*cos(x[0]-x[2])*pow(x[3], 2))
		 )/( p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2) );
    }
  };
  //]


  double get_s1(double th1, double th2) {
    static auto s1 = get_csv( "sth1dot.csv" );
    static double dth = 0.10755;
    static double min_th = -3.38782;
    AngleWrap( th1 ); AngleWrap( th2 );
    // Alt 1:
    size_t xindx = round( (th2 - min_th)/dth );
    size_t yindx = round( (th1 - min_th)/dth );
    return s1[xindx][yindx];    
    //
    // Alt 2:
    // static double frac2, frac1, intpart2, intpart1;
    // static size_t xindx_m, xindx_p, yindx_m, yindx_p;
    // frac2 = modf((th2 - min_th)/dth, &intpart2);
    // frac1 = modf((th1 - min_th)/dth, &intpart1);
    // xindx_m=intpart2; xindx_p=xindx_m+1;
    // yindx_m=intpart1; yindx_p=yindx_m+1;
    // if (xindx_p > 64 || xindx_m < 0 || yindx_p > 64 || yindx_m < 0 )
    //   return s1[round( (th2 - min_th)/dth )][round( (th1 - min_th)/dth )];
    // else {
    //   return s1[xindx_m][yindx_m]*(1-frac2)*(1-frac1)
    // 	+s1[xindx_p][yindx_m]*frac2*(1-frac1)
    // 	+s1[xindx_m][yindx_p]*(1-frac2)*frac1
    // 	+s1[xindx_p][yindx_p]*frac2*frac1;
    // }
  }

  double get_s2(double th1, double th2) {
    static auto s2 = get_csv( "sth2dot.csv" );
    static double dth = 0.10755;
    static double min_th = -3.38782;
    AngleWrap( th1 ); AngleWrap( th2 );
    // Alt 1:
    size_t xindx = round( (th2 - min_th)/dth );
    size_t yindx = round( (th1 - min_th)/dth );
    return s2[xindx][yindx];    
    //
    // Alt 2:
    // static double frac2, frac1, intpart2, intpart1;
    // static size_t xindx_m, xindx_p, yindx_m, yindx_p;
    // frac2 = modf((th2 - min_th)/dth, &intpart2);
    // frac1 = modf((th1 - min_th)/dth, &intpart1);
    // xindx_m=intpart2; xindx_p=xindx_m+1;
    // yindx_m=intpart1; yindx_p=yindx_m+1;
    // if (xindx_p > 64 || xindx_m < 0 || yindx_p > 64 || yindx_m < 0 )
    //   return s2[round( (th2 - min_th)/dth )][round( (th1 - min_th)/dth )];
    // else {
    //   return s2[xindx_m][yindx_m]*(1-frac2)*(1-frac1)
    // 	+s2[xindx_p][yindx_m]*frac2*(1-frac1)
    // 	+s2[xindx_m][yindx_p]*(1-frac2)*frac1
    // 	+s2[xindx_p][yindx_p]*frac2*frac1;
    // }
     
  }

  double s1_intp(double th1, double th2) {    
    return get_s1(th1,th2);
  }

  double s2_intp(double th1, double th2) {    
    return get_s2(th1,th2);
  }
  

// double s1_intp(double th1, double th2) {
//   static auto s1 = get_csv( "sth1dot.csv" );
//   static auto s2 = get_csv( "sth2dot.csv" );
//   static double dth = 0.10755;
//   static double min_th = -3.38782;
//   AngleWrap( th1 ); AngleWrap( th2 );
//   size_t xindx = round( (th2 - min_th)/dth );
//   size_t yindx = round( (th1 - min_th)/dth );
      
//   return s1[xindx][yindx];
// }

// double s2_intp(double th1, double th2) {
//   static auto s1 = get_csv( "sth1dot.csv" );
//   static auto s2 = get_csv( "sth2dot.csv" );
//   static double dth = 0.10755;
//   static double min_th = -3.38782;
//   AngleWrap( th1 ); AngleWrap( th2 );
//   size_t xindx = round( (th2 - min_th)/dth );
//   size_t yindx = round( (th1 - min_th)/dth );
      
//   return s2[xindx][yindx];
// }

/* 
   Note:
   dxdesdx = [ 0        ds1dth1  0        ds2dth1 ]
             [ 0        0        0        0       ]
	     [ 0        ds1dth2  0        ds2dth2 ]
	     [ 0        0        0        0       ]
*/
mat_type Gxdes(double th1, double th2) {
  static double dth = 0.108;
  static mat_type gxdesdx = mat_type::Zero(4,4);
  // gsdth1
  gxdesdx(0,1) = (s1_intp(th1+dth,th2)-s1_intp(th1-dth,th2))/(2.0*dth);
  gxdesdx(0,3) = (s2_intp(th1+dth,th2)-s2_intp(th1-dth,th2))/(2.0*dth);
  // gsdth2
  gxdesdx(2,1) = (s1_intp(th1,th2+dth)-s1_intp(th1,th2-dth))/(2.0*dth);
  gxdesdx(2,3) = (s2_intp(th1,th2+dth)-s2_intp(th1,th2-dth))/(2.0*dth);
  
  return gxdesdx;
}


/*!
    Read in a numeric csv to a 2D vector of doubles
    \param[in] filename The csv filename to read in
*/
std::vector< std::vector< double > > get_csv( const char * filename ) {
  using namespace std; 
  ifstream file( filename );
  string line;
  vector< vector< double > > mat;

  while ( getline( file, line ) ) {
    stringstream strstr( line );
    string svalue;
    vector< double > vals;
    while ( getline( strstr, svalue, ',' ) ) {
      vals.push_back( stod(svalue) );
    }
    mat.push_back( std::move(vals) );
  }
  return mat;
}

}

#endif  // SYS_DYNAM_HPP
