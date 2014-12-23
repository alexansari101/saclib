#ifndef SYS_LIN_HPP
#define SYS_LIN_HPP

namespace sac {

  //[ Linearizations of the system defined as a class
  // USER SPECIFIED:
  class sys_lin {
    double g_;
    double h_;
  
  public:
    sys_lin( ) :  g_(9.81), h_(2) {  }

    void A( const state_type & x, const state_type & u, 
	    Eigen::MatrixXd & Amat ) {
      Amat << 0, 1, 0, 0,
	( g_*cos(x[0]) - u[0]*sin(x[0]) ) / h_ ,  0, 0, 0,
	0, 0, 0, 1,
	0, 0, 0, 0;
    }

    void B( const state_type & x, const state_type & /*u*/, 
	    Eigen::MatrixXd & Bmat ) {
      Bmat << 0, cos(x[0]) / h_, 0, 1;
    }
  };
  //]

}

#endif  // SYS_LIN_HPP
