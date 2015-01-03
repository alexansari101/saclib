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
    vec_type mx_, mrho_, mrhodot_;
    mat_type mdldx_;
    mat_type mdfdx_;
  
  public:
    adjoint( state_intp & x_intp, cost & J, 
	     Params & p ) 
      :  rx_intp_( x_intp ), g_(9.81), x_(p.xlen()),
	 rho_(p.xlen()), u1_(p.ulen()), lofx_(J.m_lofx),
	 mx_(p.xlen(),1), mrho_(p.xlen(),1), mrhodot_(p.xlen(),1),
	 mdldx_(1,p.xlen()), mdfdx_(p.xlen(),p.xlen()) {  
      for ( size_t i=0; i<p.ulen(); i++ ) { u1_[i] = 0.0; } 
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
      for (indx_ = 0; indx_ < rho.size(); indx_++ ) { rhodot[indx_] = mrhodot_(indx_); }
    }
  };
  //]

}

#endif  // ADJOINT_HPP
