#ifndef SAC_STEP_HPP
#define SAC_STEP_HPP

namespace sac {

  //[ The class to carry out a sac control step
  class sac_step {
  protected:
    /* protected initialization */
    double alpha, min_val, dtt_win, dt_win;
    state_type x, rho;  state_type u_switch;  
    sys_dynam xdot;  adjoint rho_dot;
    Eigen::Matrix< double, xlen, 1 > m_mrho_tf;
    u2_optimal u2Opt;  mode_insert_grad dJdlam;
    double dJdlam_curr;   u2_cost cntrlCost;
    std::vector<double> lclMin;
    cost J1;  iter_1d it1_1d;   size_t j, steps;
  
  public:
    /* public initialization */
    double J0, Jn, t_i, t_f, t_app, tf;
    state_intp x_intp, rho_intp;  // linear interpolator
    state_type x0noU;   b_control u;
    std::vector<state_type> x_vec, rho_vec;    
    std::vector<double> times, rho_times;
    bool u2Search;  size_t its, max_its;
    Eigen::Matrix< double, xlen, 1 > m_mxdes_tf;
    
    sac_step( const bool usearch, 
	      void (*xdesFnptr) ( const double t, 
				  Eigen::Matrix< double, xlen, 1 > &m_mxdes )
	     ) : dt_win(maxdt/2.0), x(xlen), rho(xlen), 
		 u_switch(ulen), xdot(u), rho_dot(x_intp, J1),
		 u2Opt(x_intp, rho_intp, alpha),
		 dJdlam( x_intp, rho_intp, u2Opt ),
		 cntrlCost( u2Opt, dJdlam ),
		 J1(x_intp, xdesFnptr, m_mxdes_tf),
		 t_i(T), t_f(T), tf(T), 
		 x_intp(x_vec , times, xlen),
		 rho_intp(rho_vec , rho_times, xlen),
		 x0noU(xlen), u(ulen), u2Search(usearch),
                 max_its(4), m_mxdes_tf( Eigen::Matrix
					 < double, xlen, 1 >
					 ::Zero(xlen,1)) { }
    
    inline virtual void SimInitXRho( double &t0, const state_type &x0, 
				     const state_type &u_old, 
				     const double &t_i_old, 
				     const double &t_f_old );
 
    inline virtual void SimNewX( double &t0, const state_type &x0, 
				 const state_type &u_old, const double &t_i_old, 
				 const double &t_f_old );

    inline void operator() ( double &t0, const state_type &x0, 
			     const state_type &u_old, const double &t_i_old, 
			     const double &t_f_old ) {
      /* Initialize final time */
      tf = t0+T;

      /* Simulate initial trajectory */      
      SimInitXRho( t0, x0, u_old, t_i_old, t_f_old );

      /* Initialize Cost, timestep, & iters before applying u2 */
      J1.update( );
      J0 = J1;             Jn = J0;
      dtt_win = dt_win;    its = 0;

      /* Set alpha based on the initial cost */
      alpha = lam*J0;

      /* u2 automatically computed from x_intp - find opt time to apply */
      if ( u2Search ) { lclMin.clear();
	MinSearch( cntrlCost, t0+calc_tm, tf, lclMin, 0.006, 1E-3 ); 
	it1_1d = get_min(lclMin.begin(), lclMin.end(), cntrlCost, min_val);
	t_app = *it1_1d; /* OR: */ } else { t_app = t0+calc_tm; }
    
      u2Opt( t_app, u_switch );    dJdlam( t_app, dJdlam_curr );
    
      t_i = (t_app-dtt_win); t_f = (t_app+dtt_win);
      if ( t_i < t0+calc_tm ) { t_i = t0+calc_tm; }   
      if ( t_f > tf ) { t_f = tf; }   if ( t_f < t_i ) { t_i = t_f; }
      //
      else if ( dJdlam_curr < 0 ) {  // use control only if dJdlam < 0
	/* Simulate X based on applying u2* at optimal time */
	SimNewX( t0, x0, u_old, t_i_old, t_f_old );
      
	/* Update Cost after applying u2* */
	J1.update();	Jn = J1;

	/* Backtrack until cost is improved by desired amount or its */
	while ( ( (Jn>J0) || (std::abs(Jn-J0)<0.01*J0) ) && (its<max_its) ) {
	  dtt_win = dtt_win/2.0;
	  t_i = (t_app-dtt_win); t_f = (t_app+dtt_win);
	  if ( t_i < t0+calc_tm ) { t_i = t0+calc_tm; }  
	  if ( t_f > tf ) { t_f=tf; } if ( t_i>=t_f ) { t_i=t_f; break; }
	  else {  
	    /* Simulate X based on applying new duration */
	    SimNewX( t0, x0, u_old, t_i_old, t_f_old );
      
	    /* Update Cost after applying u2* */
	    J1.update();      Jn = J1;	  
	  } /* end else */
	  its++;
	} /* end backtracking while */
      
      } /* end else if */  else { Jn = J0+1; }

      if ( t_f > t0+ts+calc_tm ) { t_f = t0+ts+calc_tm; }  
      if ( t_f < t_i ) { t_i = t_f; }
    }

  };


  inline void sac_step::SimInitXRho( double &t0, const state_type &x0, 
				     const state_type &u_old, 
				     const double &t_i_old, 
				     const double &t_f_old ) {
    u=u_old; u.stimes( t_i_old, t_f_old );
    x_vec.clear(); times.clear(); x = x0;  // re-initialize state 
    simX( xdot, x, t0, t0+calc_tm, x_vec, times );  //  old_u2
    x_vec.pop_back(); times.pop_back();
    simX( xdot, x, t0+calc_tm, tf, x_vec, times );  //  no control
    x_intp( t0+ts, x0noU );
    rho_vec.clear(); rho_times.clear();  // empty the vectors
    m_mrho_tf = J1.get_dmdx( );
    for ( j = 0; j < xlen; j++ ) { rho[j] = m_mrho_tf[j]; }
    steps = simRho( rho_dot, rho, t0, tf, rho_vec, rho_times );    
  }

  inline void sac_step::SimNewX( double &t0, const state_type &x0, 
				const state_type &u_old, 
				const double &t_i_old, 
				const double &t_f_old ) {
    x_vec.clear(); times.clear(); x = x0;  // re-initialize state 
    u=u_old; u.stimes( t_i_old, t_f_old );
    simX( xdot, x, t0, t0+calc_tm, x_vec, times ); // u_old
    x_vec.pop_back(); times.pop_back();
    u=u_switch; u.stimes( t_i, t_f );              // update u
    simX( xdot, x, t0+calc_tm, t_i, x_vec, times ); // pre-u_new 
    x_vec.pop_back(); times.pop_back();
    simX( xdot, x, t_i, t_f, x_vec, times ); // u_new 
    x_vec.pop_back(); times.pop_back();
    simX( xdot, x, t_f, tf, x_vec, times ); // post-u_new 
  }
  //]  

}

#endif  // SAC_STEP_HPP
