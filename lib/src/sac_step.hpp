#ifndef SAC_STEP_HPP
#define SAC_STEP_HPP

namespace sac {

  //! The class to carry out a sac control step
  class sac_step {
  protected:
    /* protected initialization */
    double alpha_, min_val_, dtt_win_, dt_win_;
    state_type x_, rho_;  state_type u_curr_;  
    sys_dynam xdot_;  adjoint rho_dot_;
    vec_type m_mrho_tf_;
    u2_optimal u2Opt_;  mode_insert_grad dJdlam_;
    double dJdlam_curr_;   u2_cost cntrlCost_;
    std::vector<double> lclMin_;
    cost J1_;  iter_1d it1_1d_;   size_t j_, steps_;
    Params & p_;

    inline virtual void SimInitXRho( const double &t0, const state_type &x0, 
				     const b_control &u1 );

    inline virtual void SimNewX( const double &t0, const state_type &x0, 
				 const b_control &u1 );
  
  public:
    /* public initialization */
    double J0, Jn, t_i, t_f, t_app, tf;
    state_intp x_intp, rho_intp;  // linear interpolator
    state_type x0noU;   b_control u;
    std::vector<state_type> x_vec, rho_vec;    
    std::vector<double> times, rho_times;
    size_t its;
    
    sac_step( Params & p ) 
      : dt_win_(p.maxdt()/2.0), x_(p.xlen()), rho_(p.xlen()), 
	u_curr_(p.ulen()), xdot_(u), rho_dot_(x_intp, J1_, p),
	m_mrho_tf_( p.xlen(), 1 ),
	u2Opt_( x_intp, rho_intp, alpha_, p ),
	dJdlam_( x_intp, rho_intp, u2Opt_, p ),
	cntrlCost_( u2Opt_, dJdlam_, p ),
	J1_(x_intp, p), p_(p),
	t_i(p.T()), t_f(p.T()), tf(p.T()), 
	x_intp(x_vec , times, p.xlen()),
	rho_intp(rho_vec , rho_times, p.xlen()),
	x0noU(p.xlen()), u(p)
    { }

    /*!
      Performs an iteration of SAC control.
      \param[in,out] t0 initial time associated with state vector input x0.  
      This get updated by one sample time to t0=t0+ts after stepper completes.
      \param[in,out] x0 initial state vector to be integrated forward after one
      iteration of SAC.  The stepper integrates x0 from time t0 to time t0+ts 
      based on SAC controls.
      \param[in,out] u1 default control value to apply from t0 to t0+calc_tm.
      The stepper computes the new SAC control to appy from t0+calc_tm to 
      t0+calc_tm+ts and returns it in u1.
     */
    inline void operator() ( double &t0, state_type &x0, 
			     b_control &u1 ) {
      /* Initialize final time */
      tf = t0+p_.T();

      /* Simulate initial trajectory */
      SimInitXRho( t0, x0, u1 );

      /* Initialize Cost, timestep, & iters before applying u2 */
      J1_.update( );
      J0 = J1_;             Jn = J0;
      dtt_win_ = dt_win_;    its = 0;

      /* Set alpha based on the initial cost */
      alpha_ = p_.lam()*J0;

      /* u2 automatically computed from x_intp - find opt time to apply */
      if ( p_.u2search() ) { lclMin_.clear();
	MinSearch( cntrlCost_, t0+p_.calc_tm(), tf, lclMin_, 0.006, 1E-3 ); 
	it1_1d_=get_min(lclMin_.begin(), lclMin_.end(), cntrlCost_, min_val_);
	t_app = *it1_1d_; /* OR: */ } else { t_app = t0+p_.calc_tm(); }
    
      u2Opt_( t_app, u_curr_ );    dJdlam_( t_app, dJdlam_curr_ );
    
      t_i = (t_app-dtt_win_); t_f = (t_app+dtt_win_);
      if ( t_i < t0+p_.calc_tm() ) { t_i = t0+p_.calc_tm(); }   
      if ( t_f > tf ) { t_f = tf; }   if ( t_f < t_i ) { t_i = t_f; }
      //
      else if ( dJdlam_curr_ < 0 ) {  // use control only if dJdlam < 0
	/* Simulate X based on applying u2* at optimal time */
	SimNewX( t0, x0, u1 );
      
	/* Update Cost after applying u2* */
	J1_.update();	Jn = J1_;

	/* Backtrack until cost is improved by desired amount or its */
	while ( ( (Jn>J0) || (std::abs(Jn-J0)<0.01*J0) ) 
		&& (its<p_.backtrack_its()) ) {
	  dtt_win_ = dtt_win_/2.0;
	  t_i = (t_app-dtt_win_); t_f = (t_app+dtt_win_);
	  if ( t_i < t0+p_.calc_tm() ) { t_i = t0+p_.calc_tm(); }  
	  if ( t_f > tf ) { t_f=tf; } if ( t_i>=t_f ) { t_i=t_f; break; }
	  else {  
	    /* Simulate X based on applying new duration */
	    SimNewX( t0, x0, u1 );
      
	    /* Update Cost after applying u2* */
	    J1_.update();      Jn = J1_;	  
	  } /* end else */
	  its++;
	} /* end backtracking while */
      
      } /* end else if */  else { Jn = J0+1; }

      if ( t_f > t0+p_.ts()+p_.calc_tm() ) { t_f = t0+p_.ts()+p_.calc_tm(); }  
      if ( t_f < t_i ) { t_i = t_f; }
      
      if ( Jn > J0 ) { // cost increased so don't apply control
      	x0 = x0noU;    // return state under default control
	t_i=t0+p_.calc_tm(); t_app=t_i; 
      	t_f=t_i+p_.ts();    // update control horizon
	u1.clear(); u1.stimes(t_i, t_f); // return default u1 over horizon
      }
      else { 
      	x_intp( t0+p_.ts(), x0 ); // return updated state
	// return new control with new horizon {t_i, t_app, t_f}
	u1=u_curr_; u1.stimes(t_i, t_f);
      }
      t0=t0+p_.ts(); // return updated time
    }

  };

  // simulate to update x_vec, rho_vec, time vecs and interpolation objects
  inline void sac_step::SimInitXRho( const double &t0, const state_type &x0, 
				     const b_control &u1 ) {
    u=u1;
    x_vec.clear(); times.clear(); x_ = x0;  // re-initialize state 
    simX( xdot_, x_, t0, t0+p_.calc_tm(), x_vec, times );  //  old_u2
    x_vec.pop_back(); times.pop_back();
    simX( xdot_, x_, t0+p_.calc_tm(), tf, x_vec, times );  //  no control
    x_intp( t0+p_.ts(), x0noU );
    rho_vec.clear(); rho_times.clear();  // empty the vectors
    m_mrho_tf_ = J1_.get_dmdx( );
    for ( j_ = 0; j_ < rho_.size(); j_++ ) { rho_[j_] = m_mrho_tf_(j_); }
    steps_ = simRho( rho_dot_, rho_, t0, tf, rho_vec, rho_times );    
  }
  
  // simulate to update x_vec, times and interpolation objects
  inline void sac_step::SimNewX( const double &t0, const state_type &x0, 
				 const b_control &u1 ) {
    x_vec.clear(); times.clear(); x_ = x0;  // re-initialize state 
    u=u1;
    simX( xdot_, x_, t0, t0+p_.calc_tm(), x_vec, times );  // u_old
    x_vec.pop_back(); times.pop_back();
    simX( xdot_, x_, t0+p_.calc_tm(), t_i, x_vec, times ); // pre-u_new 
    x_vec.pop_back(); times.pop_back();
    u=u_curr_; u.stimes( t_i, t_f );                // update u
    simX( xdot_, x_, t_i, t_f, x_vec, times );       // u_new 
    x_vec.pop_back(); times.pop_back();
    simX( xdot_, x_, t_f, tf, x_vec, times );        // post-u_new 
  }  

}

#endif  // SAC_STEP_HPP
