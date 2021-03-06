/**
 * A header-only library implementation of sequential action control (SAC)
 *  
 * Copyright (C) 2016 Alex R. Ansari and Todd D. Murphey
 *
 * This file is part of SAClib.
 *
 * SAClib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAClib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAClib.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef TRAJ_COST_HPP
#define TRAJ_COST_HPP

namespace sac {

  /*!
    Computes the value of incremental trajectory tracking cost, \f$l(x(t))\f$ when
    called by a trajectory cost object.
  */
  class inc_state_cost {
    state_intp & rx_intp_;
    mat_type & rmQ_;
    state_type x_;
    vec_type mx_, mxdes_;
    Params & p_;

  public:

    /*!
      Constructs a inc_state_cost object to computes the value of incremental trajectory 
      tracking cost, \f$l(x(t))\f$ when called by a trajectory cost object.
      \param[in] rx_intp Reference to state interpolation object
      \param[in] rmQ  Weight matrix defining norm on incremental state 
      tracking error
      \param[in] p SAC parameters
    */
    inc_state_cost( state_intp & rx_intp,
		    mat_type & rmQ,
		    Params & p
		    ) : rx_intp_( rx_intp ),
			rmQ_( rmQ ),
			x_(p.xlen()),
			mx_( vec_type::Zero(p.xlen(),1) ),
			mxdes_( vec_type::Zero(p.xlen(),1) ),
			p_(p) { }


    /*!
      Computes the value of incremental trajectory tracking cost, \f$l(x(t))\f$.
      \param[in] J The previous cost
      \param[out] dJdt The current value of the incremental cost
      \param[in] t The current time
    */
    void operator () (const state_type &/*J*/, state_type &dJdt, const double t)
    {
      rx_intp_(t, x_); // store the current state in x
      p_.proj( x_ ); // project state if necessary
      State2Mat( x_, mx_ ); // convert state to matrix form
      //
      p_.x_des( t, x_, mxdes_ ); // Get desired trajectory point
      //
      dJdt[0] = ( ( (mx_- mxdes_).transpose() * rmQ_ 
		    * (mx_- mxdes_) ) / 2.0 )(0);
    }


  };
  //]




  /*!
    Evaluates the trajectory cost over a defined horizon from state and
    control matrices.  
    \f$J_1 = \frac{1}{2}\int_{t_0}^{t_f} (x-x_{des})^T Q\, (x-x_{des}) + u^T R\, u  \, dt + \frac{1}{2}(x-x_{des})^T P\, (x-x_{des})\f$
    Keeps references to state interpolator so that changes in state trajectory 
    are automatically accounted for.
  */
  class traj_cost {
    state_intp & rx_intp_;
    std::vector<state_type> & ru2list_, & rTiTappTf_;
    vec_type mxdes_tf_;
    mat_type mQ_;
    mat_type mP1_;
    mat_type mR_;
    state_type x_, cost_;
    vec_type mx_;
    vec_type mu2_;
    double integ_state_cost_, u2cost_, term_cost_;
    inc_state_cost inc_tracking_cost_;
    Params & p_;

  public:

    //! \todo Alex: Make constructor reference and pointer inputs const.
    /*!
      Constructs a traj_cost object from a state interpolation object, desired
      state trajectory, and the list of controls and application times.
      \param[in] rx_intp Reference to state interpolation object
      \param[in] ru2list Reference to list of applied control
      \param[in] rTiTappTf Reference to list of application times for u2
      \param[in] rmxdes_tf Reference to \f$x_{des}(t_f)\f$
      \param[in] p SAC parameters
    */
    traj_cost( state_intp & rx_intp,
	       std::vector<state_type> & ru2list,
	       std::vector<state_type> & rTiTappTf,
	       const vec_type & rmxdes_tf,
	       Params & p ) 
      : rx_intp_( rx_intp ),
	ru2list_( ru2list ),
	rTiTappTf_( rTiTappTf ),
	mxdes_tf_( rmxdes_tf ),
	mQ_( mat_type::Identity(p.xlen(),p.xlen()) ),
        mP1_( mat_type::Zero(p.xlen(),p.xlen()) ),
        mR_( mat_type::Identity(p.ulen(),p.ulen()) ),
        x_(p.xlen()),
        cost_(1),
        mx_( vec_type::Zero(p.xlen(),1) ),
        mu2_( vec_type::Zero(p.ulen(),1) ),
        integ_state_cost_(0.0), 
        u2cost_(0.0),
        term_cost_(0.0),
        inc_tracking_cost_( rx_intp, mQ_, p ),
        p_(p)  { }

    //! \todo Alex: Make constructor reference and pointer inputs const.
    /*!
      Constructs a traj_cost object from a state interpolation object, desired
      state trajectory, and the list of controls and application times.  This 
      alternate constructor assumes no terminal cost on trajectory.
      \param[in] rx_intp Reference to state interpolation object
      \param[in] ru2list Reference to list of applied control
      \param[in] rTiTappTf Reference to list of application times for u2
      \param[in] p SAC parameters
    */
    traj_cost( state_intp & rx_intp,
    	       std::vector<state_type> & ru2list,
    	       std::vector<state_type> & rTiTappTf,
    	       Params & p ) 
      : traj_cost( rx_intp, ru2list, rTiTappTf, 
    		   vec_type::Zero(p.xlen(),1), p )  { }

  
    /*!
      The function computes the integral of \f$l(x(t),u(t))\f$ and appends it to the
      terminal cost to return cost \f$J = \int_{t_0}^{t_f} l(x(t),u(t)) \,
      dt + m(x(t_f))\f$.
    */
    size_t compute_cost( const double t0, const double tf );

    /*!
      Prints the integrated control and state tracking costs, the terminal cost, and
      the total cost of the trajectory as of the last call to compute_cost( ).  Results
      are outputed to std out.
    */
    void print( ) {
      using namespace std;
      cout << "Integrated control cost: " << u2cost_ << endl;
      cout << "Integrated state tracking cost: " << integ_state_cost_ << endl;
      cout << "Terminal state cost: " << term_cost_ << endl;
      cout << "Total trajectory cost: " << cost_[0] << endl;
    }

    /*!
      \return The trajectory cost computed from the last call to compute_cost()
    */
    double get_cost( ) { return cost_[0]; }

    /*!
      \param[out] cost_vec Vector of 1) total trajectory cost 2) state cost 3) control cost 
      and 4) terminal cost
    */
    void get_costs( std::vector<double> & cost_vec ) {
      cost_vec.clear();
      cost_vec.push_back(cost_[0]);
      cost_vec.push_back(integ_state_cost_);
      cost_vec.push_back(u2cost_);
      cost_vec.push_back(term_cost_);
    }

    /*!
      get using:      mat_type rQ = J_traj.Q();
      get ref using:  mat_type & rQ = J_traj.Q();
      set using:      J_traj.Q() << 1000, 0, 0, 10;
      \return A reference to the mQ_ weight matrix for both getting and setting
    */
    mat_type & Q( ) { return mQ_; }

    /*!
      \return A reference to the mP1_ weight matrix for both getting and setting
    */
    mat_type & P( ) { return mP1_; }

    /*!
      \return A reference to the mR_ weight matrix for both getting and setting
    */
    mat_type & R( ) { return mR_; }

  };
  //]


  /*!
    \param[in] t0 The initial time for integration.
    \param[in] tf The final time for integration.
    \return The number of integration steps required.
  */
  size_t traj_cost::compute_cost( const double t0, const double tf ) {
    using namespace std;
    using namespace boost::numeric::odeint;
    double eps = 1E-7;  // if tf is == term time the floating precision can cause seg faults

    // compute terminal cost
    term_cost_ = 0;  
    rx_intp_(tf-eps, x_);
    p_.proj( x_ ); // project state if necessary
    State2Mat( x_, mx_ ); // convert state to matrix form
    cost_[0] = ( ( (mx_- mxdes_tf_).transpose() * mP1_ 
		   * (mx_- mxdes_tf_) ) / 2.0 )(0);
    term_cost_ = cost_[0];
  
    // integrate appending incremental state tracking cost to terminal cost
    typedef runge_kutta_dopri5< state_type > stepper_type;

    size_t J_steps = integrate_adaptive( make_controlled( p_.eps_cost( ) , 
							  p_.eps_cost( ) , 
							  stepper_type( ) ) , 
					 inc_tracking_cost_ , cost_ , t0 , tf-eps , 0.01 );

    integ_state_cost_ = cost_[0]-term_cost_; 

    // integrate appending incremental control cost to terminal and state tracking cost
    size_t len;
    try { len = ru2list_.size();
      if ( len != rTiTappTf_.size() ) { 
	std::cout << "Error, size of u2list does not match the number of application times in TiTappTf.\n";
	return J_steps;
      } 
    }
    catch( const std::exception& e ) { 
      std::cout << "Exception in computing the size() of ru2list_ or rTiTappTf_ in traj_cost.hpp  "
		<< e.what() <<"\n";
    }
    
    double dt; eps = 0.0000001; 
    u2cost_ = 0;
    for ( size_t i = 0; i < ru2list_.size(); i++ ) {
      dt = rTiTappTf_[i][2]-rTiTappTf_[i][0];
      if ( dt > eps) {
	State2Mat( ru2list_[i], mu2_ ); // convert state to matrix form
	u2cost_ += ( ( mu2_.transpose() * mR_ 
		       * mu2_ ) * dt )(0);
      }
    }
  
    u2cost_ = u2cost_/2.0;
    cost_[0] += u2cost_;

    return J_steps;
  }
  //]

}

#endif  // TRAJ_COST_HPP
