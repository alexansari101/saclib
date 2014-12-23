// #include "CartPend.hpp"
#ifndef CARTPEND_HPP
#define CARTPEND_HPP

#include <iostream>
#include "setting.hpp"              // Project globals & settings
#include <master.hpp>               // Master include file

using namespace sac;

/* Required by lib fuctions */
void initialize();
std::vector<double> sac_stepper(std::vector<double> xinit, double tinit);
inline void state_proj( state_type & x );  // required by traj_cost class
inline void get_DesTraj( const double t, const state_type &x,
			 Eigen::MatrixXd &m_mxdes );

namespace sac {
  namespace init {
    using namespace std;
    //[ initialization    
    state_type x(xlen), u_switch(ulen), t_curr(3);
    //
    sac_step SACit( u2Search, get_DesTraj ); 
    vector<state_type> &x_vec = SACit.x_vec, x_out, u2list, TiTappTf;
    vector<double> &times = SACit.times;
    state_intp x_intp(x_vec, times, xlen);
    state_type &x0noU = SACit.x0noU;   b_control &u = SACit.u;
    double &J0 = SACit.J0, &Jn = SACit.Jn, &t_app = SACit.t_app, 
      &t_i = SACit.t_i, &t_f = SACit.t_f;
  }
}

/* SAC STEPPER
  input:  initial state, control to apply during calculation time horizon 
          and initial time 
  return: vector [x1_new, x2_new, t_new, u_new, t_1, t_2]
  
  x1_new - intregrated theta component of cart pendulum at time t_new
  x2_new - intregrated theta_dot component of cart pendulum at time t_new
  x3_new - intregrated cart position of cart pendulum at time t_new
  x4_new - intregrated cart verlocity of cart pendulum at time t_new
  t_new - updated time = t0 + ts  (ts specified in settings)
  u_new - control applied from [t_1, t_2] which is a subset of [t0, t0+ts].
          If [t_1, t_2] is not equal to [t0, t0+ts] then the default control
          u_new=0 is applied over the remaining interval. 
  t_1 - initial time for application of the control.  t0 <= t_1 <= t0+ts
  t_2 - final time for control application.  t0 <= t_2 <= t0+ts

  WARNING: u_new is only applied when t_2-t_1 > 0, otherwise u_new=0.
  WARNING: If [t_1, t_2] is not equal to [t0, t0+ts] then u_new=0 is applied 
           over the remaining interval.
  NOTE: for speed return and input types should be changed and passed as
        references / pointers
*/
std::vector<double> sac_stepper(std::vector<double> & xinit, 
				const std::vector<double> & u_default, 
				double tinit) {
    using namespace std;
    using namespace sac::init;

    /* Perform SAC iteration - updates: J0, Jn, u, x_intp */
    SACit( tinit, xinit, u_default, tinit, tinit+calc_tm );	  	

    /* Get new u & switching times */
    u( t_i, u_switch );   
 
    if ( Jn > J0 ) { // cost increased so don't apply control
      xinit = x0noU;    // return x
      u_switch[0]=0; // return u
      t_curr[0]=tinit+calc_tm; t_curr[1]=t_curr[0]; 
      t_curr[2]=t_curr[0]+ts; // return time horizon
    }
    else { 
      x_intp( tinit+ts, xinit );
      t_curr[0]=t_i; t_curr[1]=t_app;  
      t_curr[2]=t_f; // return time horizon
    }
    
    std::vector<double> rvec; // return vec
    rvec.push_back( xinit[0] );
    rvec.push_back( xinit[1] );
    rvec.push_back( xinit[2] );
    rvec.push_back( xinit[3] );
    rvec.push_back( tinit+ts );
    rvec.push_back( u_switch[0] );
    rvec.push_back( t_curr[0] );
    rvec.push_back( t_curr[2] );
    
    return rvec;
}


/*******************************************
   Initialize */
void initialize() {
  using namespace sac::init;
  Q << 0, 0, 0, 0,  // 200            // Q weight matrix
    0, 0, 0, 0,       // 0
    0, 0, 0, 0,     // 100
    0, 0, 0, 0;      // 50
  P << 500, 0, 0, 0,    // 0              // P weight matrix
    0, 0, 0, 0,
    0, 0, 0, 0,
    0, 0, 0, 0;
  R << 0.3;           // 0.3             // R weight matrix

  // u_default[0] = 0.0;
}


/*******************************************
   Projection for calculations involving x(t) */
inline void state_proj( state_type & x ) {
  AngleWrap( x[0] );
}


/*******************************************
   outputs a point in the desired trajectory at time t */
inline void get_DesTraj( const double /*t*/,
			 const state_type &/*x*/,
			 Eigen::MatrixXd &m_mxdes ) {
  m_mxdes << 0, 0, 0, 0;
}

#endif // CARTPEND_HPP
