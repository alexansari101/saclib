#ifndef SYS_DYNAM_HPP
#define SYS_DYNAM_HPP

namespace sac {

  /* FOR NEGATIVELY SLOPED FLOOR EXAMPLES */
  // double s = -1.0/5.0; // negative floor slope

  // //[ guard equations
  // double Phiq1q1( const state_type &x ) { return x[1] - s*x[0]; }

  // mat_type DPhiq1q1( const state_type &x ) { 
  //   mat_type DPhi = mat_type::Zero(1,4);
  //   DPhi(0,0) = -s;
  //   DPhi(0,1) = 1.0;
  //   return DPhi; 
  // }
  // //]

  // //[ reset maps
  // state_type Omegaq1q1( const state_type &x ) { 
  //   state_type Omega = x;
  //   Omega[2] = ( x[2] + s*x[3] + pow(pow(s,2.0)*pow(-s*x[2]+x[3],2.0),0.5) )
  //     /( 1 + pow(s,2.0) );
  //   Omega[3] = ( pow(s,2.0)*x[2] + pow(s,3.0)*x[3] - 
  // 		 pow( pow(s,2.0)*pow(-s*x[2]+x[3],2.0) ,0.5) )
  //     /( s + pow(s,3.0) );
  //   return Omega; 
  // }
  
  // mat_type DOmegaq1q1( const state_type &x ) { 
  //   mat_type DOmega = mat_type::Identity(4,4);
  //   double c1 = (s*x[2]-x[3]);
  //   double sqrt = pow( pow(s,2.0)*pow(-c1,2.0) ,0.5);
  //   double denom1 = 1+pow(s,2.0);
  //   double denom2 = s+pow(s,3.0);
  //   //
  //   DOmega(2,2) = ( 1.0 + (pow(s,3.0)*c1)/sqrt ) / denom1;
  //   DOmega(2,3) = ( s + sqrt/(-c1) )/ denom1;
  //   //
  //   DOmega(3,2) = ( pow(s,2.0) + (pow(s,3.0)*(-c1))/sqrt ) / denom2;
  //   DOmega(3,3) = ( pow(s,3.0) + sqrt/c1 )/ denom2;
  //   //
  //   return DOmega; 
  // }

  /* FOR SINUSOIDAL FLOOR EXAMPLE */
  // //[ guard equations
  // double Phiq1q1( const state_type &x ) { return x[1] - 0.3*cos(8*x[0]); }

  // mat_type DPhiq1q1( const state_type &x ) { 
  //   mat_type DPhi = mat_type::Zero(1,4);
  //   DPhi(0,0) = 2.4*sin(8.0*x[0]);
  //   DPhi(0,1) = 1.0;
  //   return DPhi; 
  // }
  // //]

  // //[ reset maps
  // state_type Omegaq1q1( const state_type &x ) { 
  //   state_type Omega = x;
  //   Omega[2] = ((47.0 - 72.0*cos(16.0*x[0]))*x[2] + 120.0*sin(8.0*x[0])*x[3])
  //     /(-97.0 + 72.0*cos(16.0*x[0]));
  //   Omega[3] = (120.0*sin(8.0*x[0])*x[2] + (-47.0 + 72.0*cos(16.0*x[0]))*x[3])
  //     /(-97.0 + 72.0*cos(16.0*x[0]));
  //   return Omega; 
  // }
  
  // mat_type DOmegaq1q1( const state_type &x ) { 
  //   mat_type DOmega = mat_type::Identity(4,4);
  //   //
  //   DOmega(2,0) = -((960.0*cos(8.0*x[0])*(120.0*sin(8.0*x[0])*x[2] +
  // 					  (-47.0 + 72.0*cos(16.0*x[0]))*x[3]))
  // 		    / pow((97.0 - 72.0*cos(16.0*x[0])),2));
  //   DOmega(2,2) = -1.0 + 50.0/(97.0 - 72.0*cos(16.0*x[0]));
  //   DOmega(2,3) = (120.0*sin(8.0*x[0]))/(-97.0 + 72.0*cos(16.0*x[0]));
  //   //
  //   DOmega(3,0) = (960.0*((11.0*cos(8.0*x[0]) - 36.0*cos(24.0*x[0]))*x[2] + 
  // 			  60.0*sin(16.0*x[0])*x[3]))
  //     /pow( (97.0 - 72.0*cos(16.0*x[0])) ,2);
  //   DOmega(3,2) = (120.0*sin(8.0*x[0]))/(-97.0 + 72.0*cos(16.0*x[0]));
  //   DOmega(3,3) = 1.0 + 50.0/(-97.0 + 72.0*cos(16.0*x[0]));
  //   //
  //   return DOmega; 
  // }
  
  /* FOR SIMPLE BOUNCING BALL EXAMPLE */
  //[ guard equations
  double Phiq1q1( const state_type &x ) { return x[1]; }

  mat_type DPhiq1q1( const state_type &/*x*/ ) { 
    mat_type DPhi = mat_type::Zero(1,4);
    DPhi(0,1) = 1;
    return DPhi; 
  }
  //]

  //[ reset maps
  state_type Omegaq1q1( const state_type &x ) { 
    state_type Omega = x;
    Omega[3] = -Omega[3];
    return Omega; 
  }
  
  mat_type DOmegaq1q1( const state_type &/*x*/ ) { 
    mat_type DOmega = mat_type::Identity(4,4);
    DOmega(3,3) = -1;
    return DOmega; 
  }
  
  mat_type Piq1q1( const state_type &x_minus, const vec_type &f_minus,
  		   const vec_type &f_plus ) { 

    mat_type Pi = mat_type::Identity(4,4);
    mat_type tmp = mat_type::Zero(1,4);
    
    tmp = DPhiq1q1(x_minus) / ( (DPhiq1q1(x_minus) * f_minus)[0] );
    Pi = DOmegaq1q1(x_minus) * ( mat_type::Identity(4,4) - f_minus * tmp )
      + f_plus * tmp;

    return Pi; 
  }
  //]

  //[
  /*! 
    Observer class that stores states and times during integration in user 
    specified containers of type std::vector< state_type >.  Throws a runtime
    exception when a guard becomes active.
  */
  struct event_push_back
  {
    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    int m_last_sgn;
    int m_curr_sgn;

    event_push_back( std::vector< state_type > &states , 
		     std::vector< double > &times ) : 
      m_states( states ) , m_times( times ) , m_last_sgn( -2 ) ,
      m_curr_sgn( -2 ) { }
  
    void operator()( const state_type &x , double t )
    {
      m_curr_sgn = sgn( Phiq1q1(x) );
      if ( m_curr_sgn != m_last_sgn && m_last_sgn != -2 ) {
	throw t; // maybe throw updated location
      }
      m_last_sgn = m_curr_sgn;

      m_states.push_back( x );
      m_times.push_back( t );
    }
  };
  /*! 
    \fn event_push_back::event_push_back( 
    std::vector< state_type > &states , 
    std::vector< double > &times ) : 
    m_states( states ) , m_times( times )

    Initializes member variables to reference user specified containers to
    hold the states and times to be updated during integration.
    \param[in,out] states The container to hold the vector of states.
    \param[in,out] times The container to hold the vector of times.
  */
  /*! 
    \fn void event_push_back::operator()( const state_type &x , 
    double t )

    Pushes states and times during integration unless the state crosses a
    guard. If a guard is crossed an exception is thrown capturing the time
    at which the first state is on the other side of the guard.
    \param[in] x The state to push to the user specified container holding 
    the vector of states.
    \param[in] t The time to push to the user specified container holding 
    the vector of times.
  */
  //]


  //[ The rhs of xdot = f(x) defined as a class
  // USER SPECIFIED:
  class sys_dynam {
    b_control & u_;
    state_type u_curr_;
    const size_t loc_; // hybrid location

  public:
    sys_dynam( b_control & u ) : u_(u) , u_curr_(2) , loc_(1) {  }

    void operator() (const state_type &x, state_type &dxdt, const double t)
    {
      u_(t, u_curr_);
      //
      dxdt[0] = x[2];
      dxdt[1] = x[3];
      dxdt[2] = u_curr_[0];
      dxdt[3] = u_curr_[1] - 9.8;
    }

  };
  //]


  //[
  //! \todo Alex: See about making inputs const type.
  /*!
    Simulates the state of a hybrid system forward in time from an initial 
    state at t0 to the final state at tf.
    \param[in] xdot The dynamics of the system.
    \param[in,out] x0 The initial state which gets integrated to become the 
    final state.
    \param[in] t0 The initial time.
    \param[in] tf The final time.
    \param[out] x_vec The vector of states resulting from integration.
    \param[out] times The vector of times resulting from integration.
    \param[out} events The vector of times at which state reset events occur.
    \return The number of integration steps.
  */
  size_t simHybX( sys_dynam& xdot, state_type& x0, double t0, double tf,
		  std::vector<state_type>& x_vec, std::vector<double>& times,
		  std::vector<double>& events ) {  
    using namespace std;
    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5< state_type > stepper_type;

    size_t steps = 0;

    try {
      steps = integrate_adaptive( 
				 make_controlled( 1E-7 , 1E-7 , 
						  stepper_type( ) ) , 
				 xdot , x0 , t0 , tf , 0.0005 , 
				 event_push_back(x_vec , times) );      
    } catch (double t_event) {
      // interpolate from state before guard 0-crossing
      static state_type dxdtin(x0.size()), dxdtout(x0.size()), xout(x0.size());
      static vec_type v_dxdtin = vec_type::Zero(x0.size(),1);
      static stepper_type stepper;
      static double dt, eps = 1E-3;
      static int sgn_x;
      
      x0 = x_vec.back();
      sgn_x = sgn( Phiq1q1(x0) );
      // get within eps of the current side of the guard
      while ( std::abs( Phiq1q1(x0) ) > eps ) {
	xdot( x0, dxdtin, times.back() );
	State2Mat( dxdtin, v_dxdtin );

	dt = -Phiq1q1( x0 ) / ( DPhiq1q1( x0 ) * v_dxdtin )[0];

	/* NOTE: it is possible for dt < 0 if the velocity at the point before
	   guard crossing was moving away from the guard in these cases 
	   manually set a positive dt. 
	*/
	if ( dt <= 0 )
	  dt = t_event - times.back();

	stepper.do_step( xdot, x0, dxdtin, times.back(),
			 /* outputs: */ xout, dxdtout, dt );

	while ( sgn(Phiq1q1( xout )) != sgn_x ) {
	  dt = 0.9*dt;
	  stepper.do_step( xdot, x0, dxdtin, times.back(),
			   /* outputs: */ xout, dxdtout, dt );
	}

	x0 = xout;
	x_vec.push_back( x0 );
	times.push_back( times.back()+dt );
      } // end while

      // apply reset map to state before guard 0-crossing
      x0 = Omegaq1q1( x_vec.back() );
      // push event pre-reset time t^-
      events.push_back( times.back() );
      // Define the reset state as the state at t^+ (t_event + 1E-7)    
      steps+=simHybX( xdot, x0, times.back()+1E-7, tf, x_vec, times, events );
    }

    return steps;
  }
  //]


}

#endif  // SYS_DYNAM_HPP
