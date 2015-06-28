#ifndef SYS_DYNAM_HPP
#define SYS_DYNAM_HPP

namespace sac {

  double L0=1, m=80, k=14000;

  // NOTE: I don't know t_lx or t_rx for leg in air so I can't use these yet
  // /* Ground z-height */
  // double zGrnd( const double& x ) { 
  //   return 0;
  // }
  // /*! spatial derivative of the ground function at a location
  //   \param[in] x the position at which to take the derivative.
  //   \param[out] Dgrnd spatial derivative of the ground function.
  // */
  // void DzGrnd( const double &x, state_type &Dgrnd ) {
  //   std::fill(Dgrnd.begin(), Dgrnd.end(), 0.0);
  // }

  /* Left leg length */
  double Ll( const state_type &x ) {
    return pow( pow(x[0]-x[4],2)+pow(x[2]/* - zGrnd(x) */,2) ,0.5);
  }
  state_type DLl( const state_type &x ) {
    state_type dldx(x.size());
    static double Llofx;
    Llofx = Ll(x);
    dldx[0] = ( x[0] /*+ DzGrnd[0]*(zGrnd(x) - x[2]) */ - x[4] ) / Llofx;
    dldx[1] =  0; /* DzGrnd[1]*(zGrnd(x) - x[2]) / Llofx; */
    dldx[2] = /* (1-DzGrnd[2])* */(x[2]/*- zGrnd(x)*/) / Llofx;
    dldx[3] = 0; /* DzGrnd[3]*(zGrnd(x) - x[2]) / Llofx; */
    dldx[4] = ( -x[0] /*+ DzGrnd[4] (zGrnd - x[2])*/ + x[4] ) / Llofx;
    dldx[5] = 0; /* DzGrnd[5]*(zGrnd(x) - x[2]) / Llofx; */
    dldx[6] = 0; /* DzGrnd[6]*(zGrnd(x) - x[2]) / Llofx; */
    dldx[7] = 0; /* DzGrnd[7]*(zGrnd(x) - x[2]) / Llofx; */
    return dldx;
  }
  /* Right leg length */
  double Lr( const state_type &x ) {
    return pow( pow(x[5]-x[0],2)+pow(x[2]/* - zGrnd(x) */,2) ,0.5);
  }
  state_type DLr( const state_type &x ) {
    state_type dldx(x.size());
    static double Lrofx;
    Lrofx = Lr(x);
    dldx[0] = ( x[0] - x[5] /*+ DzGrnd[0]*(zGrnd(x) - x[2]) */ ) / Lrofx;
    dldx[1] =  0; /* DzGrnd[1]*(zGrnd(x) - x[2]) / Lrofx; */
    dldx[2] = /* (1-DzGrnd[2])* */(x[2]/*- zGrnd(x)*/) / Lrofx;
    dldx[3] = 0; /* DzGrnd[3]*(zGrnd(x) - x[2]) / Lrofx; */
    dldx[4] = 0; /* DzGrnd[4]*(zGrnd(x) - x[2]) / Lrofx; */
    dldx[5] = ( x[5] - x[0] /*+ DzGrnd[5] (zGrnd(x) - x[2])*/) / Lrofx;
    dldx[6] = 0; /* DzGrnd[6]*(zGrnd(x) - x[2]) / Lrofx; */
    dldx[7] = 0; /* DzGrnd[7]*(zGrnd(x) - x[2]) / Lrofx; */
    return dldx;
  }
  /* Left leg constant */
  double Cl( const state_type &x ) { return k*( L0/Ll(x) - 1); }
  state_type DCl( const state_type &x ) {
    double mult = -k*L0/pow(Ll(x),2);
    state_type dCldx = DLl(x);
    for ( size_t i=0; i<dCldx.size(); i++ )
      dCldx[i] = mult*dCldx[i];
    return dCldx;
  }
  /* Right leg constant */
  double Cr( const state_type &x ) { return k*( L0/Lr(x) - 1); }
  state_type DCr( const state_type &x ) {
    double mult = -k*L0/pow(Lr(x),2);
    state_type dCrdx = DLr(x);
    for ( size_t i=0; i<dCrdx.size(); i++ )
      dCrdx[i] = mult*dCrdx[i];
    return dCrdx;
  }

  //[ guard equations
  /* right touch down: left single support to double support */
  double Phiq1q0( const state_type &x ){return x[2]/*-zGrnd(x)*/-L0*sin(x[6]);}
  /* left take off: double support to right single support */
  double Phiq2q1( const state_type &x ) { return Ll(x) - L0; }
  /* left touch down: right single support to double support */
  double Phiq1q2( const state_type &x ){return x[2]/*-zGrnd(x)*/-L0*sin(x[6]);}
  /* right take off: double support to left single support */
  double Phiq0q1( const state_type &x ) { return Lr(x) - L0; }
  //
  /* Return a vector of the value of each guard in the current location */
  state_type Phi(const state_type &x ) { 
    state_type phiq(1);
    static int loc;
    loc = (int) (x[7]+0.5);
    //
    switch ( loc ) {
    case 0:
      phiq[0] = Phiq1q0( x );
      return phiq;
    case 1:
      phiq[0] = Phiq2q1( x );
      phiq.push_back( Phiq0q1( x ) );
      return phiq;
    case 2:
      phiq[0] = Phiq1q2( x );
      return phiq;
    default:
      std::cout << "Location #" << loc << " has not been defined\n";
      return phiq;
    }
    //
  }
  /* Return the value of a guard that transitions between locations */
  double Phi(const state_type &x, double loc_tplus ) {
    static int loc;
    static int loc_plus;
    loc = (int) (x[7]+0.5);
    loc_plus = (int) (loc_tplus+0.5);
    //
    switch ( loc ) {
    case 0:
      return Phiq1q0( x );
    case 1:
      if (loc_plus == 0)
	return Phiq0q1( x );
      else if ( loc_plus == 2 )
	return Phiq2q1( x );
      else { std::cout << "Location #" << loc_plus 
		       << " has not been defined for location at t^+\n";
	return -1; }
    case 2:
      return Phiq1q2( x );
    default:
      std::cout << "Location #" << loc << " has not been defined\n";
      return -1;
    }
    //
  }
  //
  // derivative of guard equations
  mat_type DPhiq1q0( const state_type &x ) { 
    mat_type dPhi = mat_type::Zero(1,8);
    dPhi(0,2) = 1;    dPhi(0,6) = -L0*cos(x[6]);
    return dPhi /*- DzGrnd(x)*/;
  }
  mat_type DPhiq2q1( const state_type &x ) { 
    mat_type dPhi = mat_type::Zero(1,8);
    static state_type dLldx(8); dLldx = DLl(x);
    for (size_t i=0; i<x.size(); i++) {
      dPhi(0,i) = dLldx[i];
    }
    return dPhi;
  }
  mat_type DPhiq1q2( const state_type &x ) { 
    mat_type dPhi = mat_type::Zero(1,8);
    dPhi(0,2) = 1;    dPhi(0,6) = -L0*cos(x[6]);
    return dPhi /*- DzGrnd(x)*/;
  }
  mat_type DPhiq0q1( const state_type &x ) { 
    mat_type dPhi = mat_type::Zero(1,8);
    static state_type dLrdx(8); dLrdx = DLr(x);
    for (size_t i=0; i<x.size(); i++) {
      dPhi(0,i) = dLrdx[i];
    }
    return dPhi;
  }
  /* Return the derivative the guard that transitions between 2 locations */
  mat_type DPhi(const state_type &x, const double loc_tplus ) {
    static int loc;
    static int loc_plus;
    loc = (int) (x[7]+0.5);
    loc_plus = (int) (loc_tplus+0.5);
    //
    switch ( loc ) {
    case 0:
      return DPhiq1q0( x );
    case 1:
      if (loc_plus == 0)
	return DPhiq0q1( x );
      else if ( loc_plus == 2 )
	return DPhiq2q1( x );
      else { std::cout << "Location #" << loc_plus 
		       << " has not been defined for location at t^+\n";
	return mat_type::Zero(1,8); }
    case 2:
      return DPhiq1q2( x );
    default:
      std::cout << "Location #" << loc << " has not been defined\n";
      return mat_type::Zero(1,8);
    }
    //
  }
  //]

  //[ state reset map
  state_type Omega( const state_type &x, const double loc_tplus ) { 
    state_type omega = x;
    omega[7] = loc_tplus;
    //
    switch ( (int) (x[7]+0.5) ) {
    case 0:
      if ( x[1] >= 0 ) {
	omega[5] = x[0] + L0*cos(x[6]); }
      else {
	omega[5] = x[0] - L0*cos(x[6]); }
      return omega;
    case 1:
      return omega;
    case 2:
      if ( x[1] >= 0 ) {
	omega[4] = x[0] + L0*cos(x[6]); }
      else {
	omega[4] = x[0] - L0*cos(x[6]); }
      return omega;
    default:
      std::cout << "Location #" << (int) (x[7]+0.5) 
		<< " has not been defined\n";
    }
    // 
  }  
  // derivative of state reset maps
  mat_type DOmega( const state_type &x, const double /*loc_tplus*/ ) { 
    mat_type dOmega = mat_type::Identity(8,8); 
    dOmega(7,7) = 0;
    //
    static int loc;
    loc = (int) (x[7]+0.5);
    //
    switch ( loc ) {
    case 0:
      if ( x[1] >= 0 ) { 
	dOmega(5,0) = 1; dOmega(5,6) = -L0*sin(x[6]); }
      else { 
	dOmega(5,0) = 1; dOmega(5,6) = L0*sin(x[6]); }
      return dOmega;
    case 1:
      return dOmega;
    case 2:
      if ( x[1] >= 0 ) {
	dOmega(4,0) = 1; dOmega(4,6) = -L0*sin(x[6]); }
      else {
	dOmega(4,0) = 1; dOmega(5,6) = L0*sin(x[6]); }
      return dOmega;
    default:
      std::cout << "Location #" << loc << " has not been defined\n";
    }
    // 
  }
  //]

  //[ adjoint reset map
  mat_type Pi( const state_type &x_minus, const vec_type &f_minus,
	       const vec_type &f_plus, const double loc_tplus ) { 

    mat_type pi = mat_type::Identity(8,8);
    /* static */ mat_type tmp = mat_type::Zero(1,8);
    
    tmp = DPhi(x_minus,loc_tplus)/( (DPhi(x_minus,loc_tplus)*f_minus)[0] );
    pi = DOmega(x_minus,loc_tplus)*(mat_type::Identity(8,8) - f_minus * tmp)
      + f_plus * tmp;

    if ( x_minus[7] == loc_tplus ) {
      std::cout << "Warning, locations at t_i^- and t_i^+ are the same.\n";
    }

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
	    switch ( (int) (x[7]+0.5) /* curr location */ ) {
	    case 0:
	      throw std::make_pair(t,1.0); // t_i^+ and location q_{i+1}=1
	      break;
	    case 1:
	      if (i == 0 )
		throw std::make_pair(t,2.0); // t_i^+ and location q_{i+1}=2
	      else if (i == 1)
		throw std::make_pair(t,0.0); // t_i^+ and location q_{i+1}=0
	      else
		std::cout << "unknown location q_{i+1}\n";
	      break;
	    case 2:
	      throw std::make_pair(t,1.0); // t_i^+ and location q_{i+1}=1
	      break;
	    default:
	      std::cout << "Location #" << (int)(x[7]+0.5) 
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
    //
    double Cl_;
    double Cr_;

  public:
    sys_dynam( b_control & u ) : u_(u), u_curr_(3), loc_(0), Cl_(0), Cr_(0) { }

    void operator() (const state_type &x, state_type &dxdt, const double t)
    {
      u_(t, u_curr_);
      loc_ = (int) (x[7]+0.5); // pos. floating point location to nearest int
      //
      switch ( loc_ ) {
      case 0:
	fl(x, dxdt, t); // Left leg single support dynamics
	break;
      case 1:
	fd(x, dxdt, t); // Double support dynamics
	break;
      case 2:
	fr(x, dxdt, t); // Right leg single support dynamics
	break;
      default:
	std::cout << "Location #" << loc_ << " has not been defined\n";
	break;
      }
      //
    }    
    // Left leg single support dynamics:
    void fl(const state_type &x, state_type &dxdt, const double /*t*/) {
      Cl_ = Cl(x);
      //
      dxdt[0] = x[1];
      dxdt[1] = Cl_*(x[0]-x[4]) / m;
      dxdt[2] = x[3];
      dxdt[3] = ( Cl_*(x[2]/*-zGrnd(x[4])*/) - m*9.81 ) / m;  
      dxdt[4] = 0;
      dxdt[5] = 0;
      dxdt[6] = u_curr_[0]; 
      dxdt[7] = 0;
    }
    // Double support dynamics:
    void fd(const state_type &x, state_type &dxdt, const double /*t*/) {
      Cl_ = Cl(x); Cr_ = Cr(x);
      //
      dxdt[0] = x[1];
      dxdt[1] = ( Cl_*(x[0]-x[4]) - Cr_*(x[5]-x[0]) )/m;
      dxdt[2] = x[3];
      dxdt[3] = ( Cl_*(x[2]/*-zGrnd(x[4])*/) + Cr_*(x[2]/*-zGrnd(x[5])*/) 
		  - m*9.81 )/m;
      dxdt[4] = 0;
      dxdt[5] = 0;
      dxdt[6] = u_curr_[0]; 
      dxdt[7] = 0;
    }
    // Right leg single support dynamics:
    void fr(const state_type &x, state_type &dxdt, const double /*t*/) {
      Cr_ = Cr(x);
      //
      dxdt[0] = x[1];
      dxdt[1] = -Cr_*(x[5]-x[0])/m;
      dxdt[2] = x[3];
      dxdt[3] = ( Cr_*(x[2]/*-zGrnd(x[5])*/) - m*9.81 ) / m;  
      dxdt[4] = 0;
      dxdt[5] = 0;
      dxdt[6] = u_curr_[0]; 
      dxdt[7] = 0;
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
    } catch (std::pair<double, double> post_event) {
      static double t_event;
      t_event = post_event.first;
      // interpolate from state before guard 0-crossing
      static state_type dxdtin(x0.size()), dxdtout(x0.size()), xout(x0.size());
      static vec_type v_dxdtin = vec_type::Zero(x0.size(),1);
      static stepper_type stepper;
      static double dt, eps = 1E-7;
      static int sgn_x;
      
      x0 = x_vec.back();
      sgn_x = sgn( Phi(x0, post_event.second ) );
      // get within eps of the current side of the guard
      while ( std::abs( Phi(x0,post_event.second) ) > eps ) {
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

	while ( sgn(Phi( xout, post_event.second )) != sgn_x ) {
	  dt = 0.9*dt;
	  stepper.do_step( xdot, x0, dxdtin, times.back(),
			   /* outputs: */ xout, dxdtout, dt );
	}

	// get xout to be on the other side of the event surface
	while ( sgn(Phi( xout, post_event.second )) == sgn_x ) {
	  x0 = xout;
	  dt += eps;
	  stepper.do_step( xdot, x0, dxdtin, times.back(),
			   /* outputs: */ xout, dxdtout, dt );
	}
	dt -= eps; // (x0,dt) are before event and xout is after

	// x0 = xout;
	x_vec.push_back( x0 );
	times.push_back( times.back()+dt );
      } // end while
      // // apply reset map to state before guard 0-crossing
      // x0 = Omega( x_vec.back(), post_event.second );
      // apply reset map to state just after guard 0-crossing
      x0 = Omega( xout, post_event.second );
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
