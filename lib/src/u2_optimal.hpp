#ifndef U2_OPTIMAL_HPP
#define U2_OPTIMAL_HPP

namespace sac {

  /*! 
    Stores the optimal switching control,\f$u_2^*(t)\f$.  Keeps references to 
    user maintained state and co-state trajectory interpolation objects so that 
    \f$u_2^*(t)\f$ automatically updates along with changes in trajectory.
    \f[u_2^*(t) = \; (\Lambda + R)^{-1} \, [\Lambda \, u_1(t) + h(x(t))^T \rho(t)
    \, \alpha_d]\f]
    where \f$\Lambda = h(x(t))^T \rho(t) \rho(t)^T h(x(t))\f$ and \f$h(x(t)) = 
    \frac{\partial f_1}{\partial u_1}\f$.  The dynamics, \f$f_1\f$, should be in
    control affine form.
  */
  class u2_optimal {
    state_intp & rx_intp_;
    state_intp & rrho_intp_;
    state_type x_curr_;
    state_type rho_curr_;
    double & alpha_;
    state_type u1_;
    Eigen::Matrix< double, ulen, 1 >  mrslt_, mu1_;
    Eigen::Matrix< double, xlen, 1 > mx_curr_, mrho_curr_;
    Eigen::Matrix< double, xlen, ulen > B_;
    sys_lin lin_;
    size_t i_;
  
  public:
    //! \todo Alex: make inputs const ref type
    /*!
      Initializes references to user maintained trajectory objects and the 
      desired rate of change of a trajectory tracking cost functional.
      \param[in] x_intp User maintained state interpolation object
      \param[in] rho_intp User maintained co-state interpolation object
      \param[in] alpha User specified desired change in cost functional
      relative to the duration of activiation of \f$u_2^*(t)\f$. 
      i.e. \f$\frac{\Delta J_1}{\Delta t}\f$.
    */
    u2_optimal( state_intp & x_intp, 
		state_intp & rho_intp, 
		double & alpha ) : rx_intp_( x_intp ) ,
				   rrho_intp_( rho_intp ) ,
				   x_curr_(xlen) , rho_curr_(xlen) , 
				   alpha_( alpha ) , u1_(ulen) { 
      for ( i_=0; i_<ulen; i_++ ) {
	u1_[i_] = 0.0;
	mu1_(i_,0) = 0.0;
	mrslt_(i_,0) = 0.0;
      }
    }
  
    /*!
      Computes the value of the optimal switching control,\f$u_2^*(t)\f$.
      \f[u_2^*(t) = \; (\Lambda + R)^{-1} \, [\Lambda \, u_1(t) + h(x(t))^T 
      \rho(t) \, \alpha_d]\f]
      \param[in] t The time at which to compute the mode insertion gradient
      \param[out] u2Opt_curr The value of the optimal switching control
    */
    void operator() ( const double t, state_type & u2Opt_curr ) {
      rx_intp_(t, x_curr_);      State2Mat( x_curr_, mx_curr_ );
      rrho_intp_(t, rho_curr_);  State2Mat( rho_curr_, mrho_curr_ );
      lin_.B( x_curr_, u1_, B_ );
      mrslt_ = B_.transpose()*mrho_curr_;
      mrslt_ = ( mrslt_*mrslt_.transpose() + R ).inverse()
	*(/* (...)*mu1_ +*/ mrslt_*alpha_);    
      //
      for ( i_=0; i_<ulen; i_++ ) { u2Opt_curr[i_] = mrslt_(i_,0); }
    }
  };

}

#endif  // U2_OPTIMAL_HPP
