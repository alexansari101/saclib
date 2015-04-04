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
    std::vector<double> times, rho_times, events;
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
    u=u1; events.clear();
    x_vec.clear(); times.clear(); x_ = x0;  // re-initialize state 
    simHybX( xdot_, x_, t0, t0+p_.calc_tm(), x_vec, times, events ); // old_u2
    x_vec.pop_back(); times.pop_back(); // no control
    simHybX( xdot_, x_, t0+p_.calc_tm(), tf, x_vec, times, events );
    x_intp( t0+p_.ts(), x0noU );
    rho_vec.clear(); rho_times.clear();  // empty the vectors
    J1_.grad_mofx( m_mrho_tf_ );
    Mat2State( m_mrho_tf_, rho_ );

    if (events.size() < 1) // if no events
      steps_ = simRho( rho_dot_, rho_, t0, tf, rho_vec, rho_times );
    else {
      events.push_back( tf ); // now there are at least 2 elements in events

      while ( events.size() > 1 ) { // while there are at least 2

	/* integrate from tf (or last t_event^-) to t_event^+ */
    	steps_ = simRho( rho_dot_, rho_, *(events.rbegin()+1)+1E-7, 
    			 events.back(), rho_vec, rho_times );

	/* reverse rho and times because they grow backwards in time */
    	reverse( rho_vec.begin(), rho_vec.end() );
    	reverse( rho_times.begin(), rho_times.end() );
	
	/* apply reset map to rho */
	static state_type s_f_minus(x0.size()), s_f_plus(x0.size()),
	  s_x_minus(x0.size()), s_x_plus(x0.size());
	static vec_type v_f_minus(4,1), v_f_plus(4,1), v_rho_plus(4,1);
	// get x(t_event^-) and x(t_event^+)
	x_intp( *(events.rbegin()+1), s_x_minus );
	x_intp( *(events.rbegin()+1)+1E-7, s_x_plus );
	// get f(t_event^-) and f(t_event^+)
	xdot_( s_x_minus, s_f_minus, *(events.rbegin()+1) );
	xdot_( s_x_plus, s_f_plus, *(events.rbegin()+1)+1E-7 );
	//
	State2Mat( s_f_minus, v_f_minus );
	State2Mat( s_f_plus, v_f_plus );
	// Apply reset map to get rho(t_event^-) from rho(t_event^+)
	State2Mat( rho_, v_rho_plus );
	v_rho_plus = (Piq1q1( s_x_minus, v_f_minus, v_f_plus )).transpose()
	  *v_rho_plus; // after this it is v_rho_minus
	// update rho_ to be v_rho_plus
	Mat2State( v_rho_plus, rho_ );

    	events.pop_back();
      } // end for

      /* last event in list */
      steps_ = simRho( rho_dot_, rho_, t0, events.back(), 
    		       rho_vec, rho_times );
      events.pop_back();
    } // end else
    //]

  }
  
  // simulate to update x_vec, times and interpolation objects
  inline void sac_step::SimNewX( const double &t0, const state_type &x0, 
				 const b_control &u1 ) {

    x_vec.clear(); times.clear(); x_ = x0;  // re-initialize state 
    u=u1; events.clear();
    simHybX( xdot_, x_, t0, t0+p_.calc_tm(), x_vec, times, events );  // u_old
    x_vec.pop_back(); times.pop_back(); // pre-u_new 
    simHybX( xdot_, x_, t0+p_.calc_tm(), t_i, x_vec, times, events );
    x_vec.pop_back(); times.pop_back();
    u=u_curr_; u.stimes( t_i, t_f );                          // update u
    simHybX( xdot_, x_, t_i, t_f, x_vec, times, events );     // u_new 
    x_vec.pop_back(); times.pop_back();
    simHybX( xdot_, x_, t_f, tf, x_vec, times, events );      // post-u_new 
  }  

}

#endif  // SAC_STEP_HPP
