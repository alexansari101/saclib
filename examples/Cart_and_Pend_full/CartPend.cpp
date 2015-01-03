/*
 Full state cart pendulum system under acceleration control
 */

/******************************************************************/
/* Reserve max matrix NxN dims if known at compile time for speed */
/* Must at least be the max(x_len,u_len).                         */
// #define CONST_MAT 4

#include <iostream>
#include <master.hpp>               // Master include file
#include <boost/timer/timer.hpp>    // timing stuff 
#include <traj_cost.hpp>            // Add to compute final trajectory cost

using namespace sac;

inline void state_proj( state_type & x );  // required by traj_cost class
inline void get_DesTraj( const double t, const state_type &x,
			 vec_type &m_mxdes );

int main(int /* argc */ , char** /* argv */ )
{
    using namespace std;

    boost::timer::auto_cpu_timer t(4, "%w seconds\n");

    /* initialize SAC parameters */
    Params params(4,1);
    params.T() = 1.2;
    params.lam() = -5;
    params.maxdt() = 0.2;
    params.ts() = 0.0167;
    params.usat() = { {10, -10} };
    params.calc_tm() = 0; // params.ts();
    params.u2search() = true;
    params.Q() = mat_type::Zero(4,4);
    params.P() = mat_type::Zero(4,4);
    params.Q()(0,0) = 200; params.Q()(2,2) = 100; params.Q()(3,3) = 50;
    params.R() << 0.3;

    /* initialize SAC stepper */
    sac_step SACit( params, get_DesTraj );

    /* simulation start/stop times */
    double t0=0.0, tsim = 10;
    
    /* initial state and control */
    state_type x0(params.xlen());    x0 = { PI, 0, 0, 0 };
    b_control u1(params);     u1.stimes( 0, params.calc_tm() );

    /* compute final trajectory cost */
    vector<state_type> x_out, u2list, TiTappTf;
    vector<double> state_times;
    state_intp traj_intp(x_out, state_times, params.xlen());
    traj_cost J_traj( traj_intp, state_proj, u2list, TiTappTf, 
		      get_DesTraj, SACit.m_mxdes_tf, params );
    //
    J_traj.Q() = params.Q(); J_traj.P() = params.P();
    J_traj.R() << 0.3;
    //]
    
    /* Save initial x and u */
    x_out.push_back(x0); state_times.push_back(t0);
    u2list.push_back( u1.m_u_switch ); 
    TiTappTf.push_back( { u1.m_tau1, 0.0, u1.m_tau2 } );

    while ( t0 < tsim ) {
      
      /* Perform SAC step */
      SACit( t0, x0, u1 ); // update: u1 and x0 to time t0+ts

      /* Save updated x and u */
      x_out.push_back(x0); state_times.push_back(t0);
      u2list.push_back( u1.m_u_switch ); 
      TiTappTf.push_back( { SACit.t_i, SACit.t_app, SACit.t_f } );

    } /* end main while() loop */

    t.stop(); t.report();

    /* Save to file */
    SaveVec( x_out, "x.csv" );   SaveVec( u2list, "u2list.csv" );
    SaveVec( TiTappTf, "TiTappTf.csv" );

    //[ /* Compute final trajectory cost */
    J_traj.compute_cost( 0.0, state_times.back() );
    J_traj.print();
    //]
  
    return 0;
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
			 vec_type &m_mxdes ) {
  m_mxdes << 0, 0, 0, 0;
}
