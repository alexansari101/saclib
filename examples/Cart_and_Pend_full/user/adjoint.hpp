#ifndef ADJOINT_HPP
#define ADJOINT_HPP

namespace sac {

  //[ The rhs of rho' defined as a class
  // USER SPECIFIED:
  class adjoint {
    state_intp & rx_intp_;
    double g_, indx_;
    state_type x_, rho_, u1_;
    sys_lin lin_;
    inc_cost & lofx_;
    Eigen::MatrixXd mx_, mrho_, mrhodot_;
    Eigen::MatrixXd mdldx_;
    Eigen::MatrixXd mdfdx_;
  
  public:
    adjoint( state_intp & x_intp,
	     cost & J ) :  rx_intp_( x_intp ), g_(9.81), x_(xlen),
			   rho_(xlen), u1_(ulen), lofx_(J.m_lofx),
			   mx_(xlen,1), mrho_(xlen,1), mrhodot_(xlen,1),
			   mdldx_(1,xlen), mdfdx_(xlen,xlen) {  
      for ( size_t i=0; i<ulen; i++ ) { u1_[i] = 0.0; } 
    }

    void operator() (const state_type &rho, state_type &rhodot, const double t)
    {
      rho_ = rho;
      rx_intp_(t, x_);        // store the current state in x
      State2Mat( x_, mx_ );   // convert state to matrix form
      State2Mat( rho_, mrho_ );
      //
      lin_.A( x_, u1_, mdfdx_ );
      //
      AngleWrap( mx_(0) ); // Only for angle wrapping
      //
      lofx_.dx( t, mx_, mdldx_ );
      //
      mrhodot_ = -mdldx_.transpose() - mdfdx_.transpose()*mrho_;
      //
      for (indx_ = 0; indx_ < xlen; indx_++ ) { rhodot[indx_] = mrhodot_(indx_); }
    }
  };
  //]

}

#endif  // ADJOINT_HPP
