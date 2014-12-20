#ifndef U2_COST_HPP
#define U2_COST_HPP

namespace sac {

  /*! 
    Provides a cost function of time that can be optimized to find the best time 
    to apply the optimal control law u2*. e.g. \f$cost(t) = \sqrt{u_2^*(t)^T 
    u_2^*(t)} + \frac{dJ_1}{d \lambda^+}(t) + t^{1.6}\f$
  */
  class u2_cost {
    u2_optimal & u2Opt_;
    mode_insert_grad & dJdlam_;
    state_type u2Opt_curr_;
    double dJdlam_curr_;
    size_t i_;
    double u2_sum_;
  
  public:
    //! \todo Alex: see about converting inputs to const ref types.
    /*! 
      Initializes references to user defined objects \f$u_2^*(t)\f$ and 
      \f$\frac{dJ_1}{d \lambda^+}(t)\f$ which automatically update along with
      updates in trajectory.
      \param[in] u2Opt The object storing the values of control \f$u_2^*(t)\f$
      \param[in] dJdlam The object storing the values of mode insertion gradient 
      \f$\frac{dJ_1}{d \lambda^+}(t)\f$
    */
    u2_cost( u2_optimal & u2Opt, 
	     mode_insert_grad & dJdlam ) : u2Opt_( u2Opt ) , 
					   dJdlam_( dJdlam ),
					   u2Opt_curr_(ulen), 
					   dJdlam_curr_(0.0) { }
  
    /*! 
      Returns the value of a cost function at the specified time. The cost 
      function can be searched to find the best time to apply \f$u_2^*(t)\f$.
      \f[cost(t) = \sqrt{u_2^*(t)^T u_2^*(t)} + \frac{dJ_1}{d \lambda^+}(t) 
      + t^{1.6}\f]
      \param[in] t The time at which to evaluate the cost
      \return The value of the cost at the specified time
    */
    double operator() ( const double t )
    {
      u2Opt_(t, u2Opt_curr_);
      u2_sum_ = 0;    
      for ( i_ = 0; i_ < ulen; i_++ ) { u2_sum_ += pow(u2Opt_curr_[i_], 2); }
      u2_sum_ = pow(u2_sum_, 0.5);
      //
      dJdlam_(t, dJdlam_curr_);
      //
      return u2_sum_ + dJdlam_curr_ + pow(t, 1.6);
    }
  };

}

#endif  // U2_COST_HPP
