#ifndef SYS_DYNAM_HPP
#define SYS_DYNAM_HPP

namespace sac {

  double  L0=1, m=1, k=100;

  // Smooth Step (logistic function)
  double slope = 75.0;
  inline double smStep( const double x, const double cntr ) {
    return 1.0 /( 1.0 + exp(-slope*(x-cntr)) );
  }
  inline double DsmStep( const double x, const double cntr ) {
    return slope*exp(-slope*(x-cntr)) / pow( 1 + exp(-slope*(x-cntr)) , 2.0);
  }
  //
  

  /*!
    computes the ground z-height at a location (x,y)
    \param[in] x position in the x-dimension of the world frame.
    \param[in] y position in the y-dimension of the world frame.
  */
  inline double zGrndToe( const double& x, const double& y ) {     
    // return .5*x+0.3*cos(4.0*x - PI)*cos(4.0*y - PI); // 0.3*cos(4.0*x - PI);
    return .5*( smStep(x,.7) + smStep(x,1.4) + smStep(x,2.1) + smStep(x,2.8) ) ;
    // return 0;
  }
  /*!
    computes the spatial derivative of the ground function at a location (x,y)
    \param[in] x state_type storing the current state.
    \param[out] Dgrnd spatial derivative of the ground function.
  */
  inline void DzGrndToe( const state_type &x, state_type &Dgrnd ) {
    for (size_t i = 0; i<9; i++) {
      Dgrnd[i]=0;
    }
    // NOTE!!! This used to be Dgrnd[4]= Fun(x[4]).. 
    //         but zGrnd is fun of x[6] and x[7] only!
    // Dgrnd[6]= /*.5*/-1.2*sin(4.0*x[6]-PI); // Dgrnd[6]= -1.2*sin(4.0*x[6]-PI);
    Dgrnd[6] = .5*( DsmStep(x[6],.7) + DsmStep(x[6],1.4) 
		    + DsmStep(x[6],2.1) + DsmStep(x[6],2.8) ) ;
    // double finite difference
    // static double dx = .1;
    // Dgrnd[6]= (zGrndToe(x[6]+dx,x[7])-zGrndToe(x[6]-dx,x[7]))/(2.0*dx);
  }

  //[ guard equations
  /* left take off: left single support to aerial */
  double Phiq3q0( const state_type &x ) { 
    // // TEST
    // state_type xx(x.size());
    // for (size_t i=0; i<x.size(); i++) {
    //   xx[i]=x[i];
    //   if (xx[i] < 1E-16){
    //   	  xx[i] = 0;
    //   	}
    //   }
    // return xx[4] // TEST
    //   - L0*( xx[4]-zGrndToe(xx[6], xx[7]) )
    //   / pow( pow(xx[0]-xx[6],2)+pow(xx[2]-xx[7],2)+pow(xx[4]-zGrndToe(xx[6],xx[7]),2) 
    // 	     ,0.5)
    //   - zGrndToe(xx[6],xx[7]);  
    //
    return x[4] 
      - L0*( x[4]-zGrndToe(x[6], x[7]) )
      / pow( pow(x[0]-x[6],2)+pow(x[2]-x[7],2)+pow(x[4]-zGrndToe(x[6],x[7]),2) 
    	     ,0.5)
      - zGrndToe(x[6],x[7]);  
  }
  /* left touch down: aerial to left single support */
  double Phiq0q3( const state_type &x ) { return Phiq3q0(x); }
  /**/
  //
  /* Return a vector of the value of each guard in the current location */
  state_type Phi(const state_type &x ) { 
    state_type phiq(1);
    static int loc;
    loc = (int) (x[8]+0.5);
    //
    switch ( loc ) {
    case 0:
      phiq[0] = Phiq3q0( x );
      return phiq;
    case 3:
      phiq[0] = Phiq0q3( x );
      return phiq;
    default:
      std::cout << "Phi: Location #" << loc << " has not been defined\n";
      return phiq;
    }
    //
  }
  /* Return the value of a guard that transitions between locations */
  double Phi(const state_type &x, double loc_tplus ) {
    static int loc;
    static int loc_plus;
    loc = (int) (x[8]+0.5);
    loc_plus = (int) (loc_tplus+0.5);
    //
    switch ( loc ) {
    case 0:
      if (loc_plus == 3)
        return Phiq3q0( x );
      else { std::cout << "Phi: Location #" << loc_plus 
		       << " has not been defined for location at t^+\n"; 
        return -1; }
    case 3:
      if (loc_plus == 0)
	return Phiq0q3( x );
      else { std::cout << "Phi: Location #" << loc_plus 
		       << " has not been defined for location at t^+\n";
      	return -1; }
    default:
      std::cout << "Phi: Location #" << loc << " has not been defined\n";
      return -1;
    }
    //
  }
  //
  mat_type DPhiq3q0( const state_type &x ) { 
    mat_type dPhi = mat_type::Zero(1,9);
    static state_type Dgrnd(9);
    DzGrndToe( x, Dgrnd );
    //
    double denom32 = pow( pow(x[0] - x[6],2.0) + pow(x[2] - x[7],2.0) 
			  + pow(x[4] - zGrndToe(x[6], x[7]),2), 3.0/2.0);
    //
    dPhi(0,0) = (L0*(x[0] - x[6])*(x[4] - zGrndToe(x[6], x[7]))) / denom32;
    dPhi(0,2) = (L0*(x[2] - x[7])*(x[4] - zGrndToe(x[6], x[7]))) / denom32;
    dPhi(0,4) = 1 - L0 / pow( pow(x[0] - x[6],2.0) + pow(x[2] - x[7],2.0) 
			      + pow(x[4] - zGrndToe(x[6], x[7]),2), 0.5)  
      + (L0*pow(x[4] - zGrndToe(x[6], x[7]),2)) / denom32;
    //
    dPhi(0,6) = -Dgrnd[6] + L0*Dgrnd[6] / 
      pow( pow(x[0] - x[6],2.0) + pow(x[2] - x[7],2.0) 
	   + pow(x[4] - zGrndToe(x[6], x[7]),2), 0.5) 
      + L0*(x[4] - zGrndToe(x[6], x[7]))
      *(-2.0*(x[0]-x[6])-2.0*(x[4] - zGrndToe(x[6], x[7]))
	*Dgrnd[6])/( 2.0*denom32 );
    dPhi(0,7) = -Dgrnd[7] + L0*Dgrnd[7] / 
      pow( pow(x[0] - x[6],2.0) + pow(x[2] - x[7],2.0) 
	   + pow(x[4] - zGrndToe(x[6], x[7]),2), 0.5) 
      + L0*(x[4] - zGrndToe(x[6], x[7]))
      *(-2.0*(x[2]-x[7])-2.0*(x[4] - zGrndToe(x[6], x[7]))
	*Dgrnd[7])/( 2.0*denom32 );
    return dPhi;
  }
  mat_type DPhiq0q3( const state_type &x ) { 
    return DPhiq3q0(x);
  }
  //]
  /**/
  /* Return the derivative the guard that transitions between 2 locations */
  mat_type DPhi(const state_type &x, const double loc_tplus ) {
    static int loc;
    static int loc_plus;
    loc = (int) (x[8]+0.5);
    loc_plus = (int) (loc_tplus+0.5);
    //
    switch ( loc ) {
    case 0:
      if (loc_plus == 3)
        return DPhiq3q0( x );
      else { std::cout << "DPhi: Location #" << loc_plus 
		       << " has not been defined for location at t^+\n";
	return mat_type::Zero(1,8); }
    case 3:
      if (loc_plus == 0)
	return DPhiq0q3( x );
      else { std::cout << "DPhi: Location #" << loc_plus 
		       << " has not been defined for location at t^+\n";
	return mat_type::Zero(1,8); }
    default:
      std::cout << "DPhi: Location #" << loc << " has not been defined\n";
      return mat_type::Zero(1,8);
    }
    //
  }
  //]

  //[ state reset map
  state_type Omega( const state_type &x, const double loc_tplus ) { 
    state_type omega = x;
    omega[8] = loc_tplus;
    //
    switch ( (int) (x[8]+0.5) ) {
    case 0:
      return omega;
    case 3:
      return omega;
    default:
      std::cout << "Omega: Location #" << (int) (x[8]+0.5) 
		<< " has not been defined\n";
      return omega;
    }
    // 
  }  
  
  // derivative of state reset maps
  mat_type DOmega( const state_type &x, const double loc_tplus ) { 
    mat_type dOmega = mat_type::Identity(9,9); 
    dOmega(8,8) = 0;
    //
    static int loc;
    loc = (int) (x[8]+0.5);
    //
    switch ( loc ) {
    case 0:
      return dOmega;
    case 3:
      return dOmega;
    default:
      std::cout << "DOmega: Location #" << loc << " has not been defined\n";
      return dOmega;
    }
    // 
  }
  //]
  
  //[ adjoint reset map
  mat_type Pi( const state_type &x_minus, const vec_type &f_minus,
	       const vec_type &f_plus, const double loc_tplus ) { 

    mat_type pi = mat_type::Identity(9,9);
    mat_type tmp = mat_type::Zero(1,9);

    if ( (int) (x_minus[8]) == (int) (loc_tplus+.05) ) {
      std::cout << "Pi: Warning, locations at t_i^- and t_i^+ are the same.\n";
    }
    
    tmp = DPhi(x_minus,loc_tplus)/( (DPhi(x_minus,loc_tplus)*f_minus)[0] );
    pi = DOmega(x_minus,loc_tplus)*(mat_type::Identity(9,9) - f_minus * tmp)
      + f_plus * tmp;

    // NOTE: f_minus and f_plus are not guaranteed / checked to be on 
    // opposite sides of the event surface.

    return pi; 
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
    std::vector< int > m_last_sgn;
    std::vector< int > m_curr_sgn;
    state_type Phi_;

    event_push_back( std::vector< state_type > &states , 
		     std::vector< double > &times ) : 
      m_states( states ) , m_times( times ) { }
  
    void operator()( const state_type &x , double t )
    {
      Phi_ = Phi(x);
      m_curr_sgn.resize( Phi_.size() );
      for (size_t i=0; i<Phi_.size(); i++)
	m_curr_sgn[i] = sgn( Phi_[i] );
      if ( m_last_sgn.size() == m_curr_sgn.size() ) {
	for ( size_t i=0; i<m_curr_sgn.size(); i++ ) {
	  if ( m_curr_sgn[i] != m_last_sgn[i] ) {
	    //
	    switch ( (int) (x[8]+0.5) /* curr location */ ) {
	    case 0:
	      // make sure we are moving away from event surface and z > 0
	      if ( x[5] > 0 && x[4] > zGrndToe(x[0],x[2]) ) { 
		throw std::make_pair(t,3.0); // t_i^+ and location q_{i+1}=3
	      }
	      break;
	    case 3:
	      // make sure we are moving away from event surface and z > 0
	      if ( x[5] < 0 && x[4] > zGrndToe(x[0],x[2]) ) {
		throw std::make_pair(t,0.0); // t_i^+ and location q_{i+1}=0
	      }
	      break;
	    default:
	      std::cout << "event_push_back: Location #" << (int)(x[8]+0.5) 
			<< " has not been defined\n";
	      break;
	    } // end switch
	  } // end inner if
	} // end for
      } // end outer if
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
    int loc_;

  public:
    sys_dynam( b_control & u ) : u_(u), u_curr_(3), loc_(0) { }

    void operator() (const state_type &x, state_type &dxdt, const double t)
    {
      u_(t, u_curr_);
      loc_ = (int) (x[8]+0.5); // pos. floating point location to nearest int
      //
      switch ( loc_ ) {
      case 0:
	fl(x, dxdt, t); // Left leg single support dynamics
	break;
      case 3:
      	ff(x, dxdt, t); // flight dynamics
      	break;
      default:
	std::cout << "sys_dynam: Location #" << loc_ << " has not been defined\n";
	break;
      }
      //
    }    
    // Left leg single support dynamics:
    void fl(const state_type &x, state_type &dxdt, const double /*t*/) {
      double L = pow( pow(x[0]-x[6],2)+pow(x[2]-x[7],2)
		      +pow(x[4]-zGrndToe(x[6],x[7]),2) ,0.5);
      dxdt[0] = x[1];
      dxdt[1] = (k/m*(L0 - L) + u_curr_[2]/m)*(x[0]-x[6])/L;
      dxdt[2] = x[3]; 
      dxdt[3] = (k/m*(L0 - L) + u_curr_[2]/m)*(x[2]-x[7])/L;
      dxdt[4] = x[5];
      dxdt[5] = (k/m*(L0 - L) + u_curr_[2]/m)*(x[4]-zGrndToe(x[6],x[7]))/L 
	- 9.81;
      dxdt[6] = 0; 
      dxdt[7] = 0;
      dxdt[8] = 0;
    }
    // flight dynamics:
    void ff(const state_type &x, state_type &dxdt, const double /*t*/) {
      dxdt[0] = x[1];
      dxdt[1] = 0;
      dxdt[2] = x[3]; 
      dxdt[3] = 0;
      dxdt[4] = x[5];
      dxdt[5] = -9.81;
      dxdt[6] = x[1]+u_curr_[0]; 
      dxdt[7] = x[3]+u_curr_[1];
      dxdt[8] = 0;
    }
  };
  //]


  //[
  //! \todo Alex: Develop a more efficient search (e.g. Golden Sect., Grad.).
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
    // std::cout << "called\n"; // TEST
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
    } catch (std::pair<double, double> post_event) {
      // std::cout << "started catch\n"; // TEST
      static double t_event;
      t_event = post_event.first;      
      // interpolate from state before guard 0-crossing
      static state_type dxdtin(x0.size()), dxdtout(x0.size()), xout(x0.size());
      static vec_type v_dxdtin = vec_type::Zero(x0.size(),1);
      static stepper_type stepper;
      static double dt, eps = 1E-6;
      static int sgn_x;
      x0 = x_vec.back();
      sgn_x = sgn( Phi(x0, post_event.second ) );
      // get within eps of the current side of the guard
      while ( std::abs( Phi(x0,post_event.second) ) > eps ) {
	// std::cout << "\n\nbeginning while.\n "; // TEST
	xdot( x0, dxdtin, times.back() );
	State2Mat( dxdtin, v_dxdtin );
	
	dt = -Phi( x0,post_event.second ) / 
	  ( DPhi( x0, post_event.second ) * v_dxdtin )[0];

	/* NOTE: it is possible for dt < 0 if the velocity at the point before
	   guard crossing was moving away from the guard in these cases 
	   manually set a positive dt. 
	*/
	if ( dt <= 0 )
	  dt = t_event - times.back();

	stepper.do_step( xdot, x0, dxdtin, times.back(),
			 /* outputs: */ xout, dxdtout, dt );

	// // TEST
	// std::cout << "x0 = (";
	// for (auto &i : x0) {
	//   cout << i << ", ";
	// }
	// std::cout << ")\n";
	// std::cout << "xout = (";
	// for (auto &i : xout) {
	//   cout << i << ", ";
	// }	  
	// std::cout << ")\n";
	// std::cout << "x differences = (";
	// for (size_t p=0; p<x0.size(); p++)
	//   std::cout << x0[p]-xout[p] << ", ";
	// std::cout << ")\n";
	// std::cout << "phi(x0) = " << Phi(x0,post_event.second)
	// 	  << "\t phi(xout) = " << Phi(xout,post_event.second) << "\n";

	// int count = 0; // TEST
	while ( sgn(Phi( xout, post_event.second )) != sgn_x ) {
	  // count++; // TEST
	  dt = 0.66*dt;
	  stepper.do_step( xdot, x0, dxdtin, times.back(),
			   /* outputs: */ xout, dxdtout, dt );
	  // for (auto &i : xout) { // TEST
	  //   if (i < 1E-8){
	  //     i = 0;
	  //   }
	  // }
	}

	// // TEST
	// if (count >= 10000) {
	//   std::cout << "\t count of 10k\n";
	//   std::cout << "x0 = (";
	//   for (auto &i : x0) {
	//     cout << i << ", ";
	//   }
	//   std::cout << ")\n";
	//   std::cout << "xout = (";
	//   for (auto &i : xout) {
	//     std::cout << i << ", ";
	//   }	  
	//   std::cout << ")\n";
	//   std::cout << "x differences = (";
	//   for (size_t p=0; p<x0.size(); p++)
	//     std::cout << x0[p]-xout[p] << ", ";
	//   std::cout << ")\n";
	//   std::cout << "phi(x0) = " << Phi(x0,post_event.second)
	// 	    << "\t phi(xout) = " << Phi(xout,post_event.second) << "\n";
	//   std::cout << ")\n\n Exiting\n\n";
	//   exit(-1);
	// }
	// std::cout << "exited ok.\n"; // TEST

	// get xout to be on the other side of the event surface
	while ( sgn(Phi( xout, post_event.second )) == sgn_x ) {
	  x0 = xout;
	  dt += eps;
	  stepper.do_step( xdot, x0, dxdtin, times.back(),
			   /* outputs: */ xout, dxdtout, dt );
	}
	dt -= eps; // (x0,dt) are before event and xout is after

	x_vec.push_back( x0 );
	times.push_back( times.back()+dt );
	// std::cout << "end while..\n"; // TEST
      } // end while
      // std::cout << "ended catch\n"; // TEST
      // apply reset map to state just after guard 0-crossing
      x0 = Omega( xout, post_event.second );
      // push event pre-reset time t^-
      events.push_back( times.back() );
      // Define the reset state as the state at t^+ (t_event + 1E-6)    
      steps+=simHybX( xdot, x0, times.back()+1E-6, tf, x_vec, times, events );
    }
    // std::cout << "returning\n"; // TEST
    return steps;
  }
  //]


}

#endif  // SYS_DYNAM_HPP
