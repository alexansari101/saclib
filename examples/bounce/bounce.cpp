/*
  Impacting point mass under acceleration control
 */

/******************************************************************/
/* Reserve max matrix NxN dims if known at compile time for speed */
/* Must at least be the max(x_len,u_len).                         */
#define CONST_MAT 4

#include <iostream>
#include <master.hpp>               // Master include file
#include <boost/timer/timer.hpp>    // timing stuff 
#include <traj_cost.hpp>            // Add to compute final trajectory cost

using namespace sac;


int main(int /* argc */ , char** /* argv */ )
{
    using namespace std;

    boost::timer::auto_cpu_timer t(4, "%w seconds\n");

    /* initialize SAC parameters */
    Params params(4,2);
    params.T() = 0.5;
    params.lam() = -10;
    params.maxdt() = 0.2;
    params.ts() = 0.01;
    params.usat() = { {10, -10}, {0, -10} };
    // params.usat() = { {0, -0}, {0, -0} };
    params.calc_tm() = params.ts();
    params.u2search() = false;
    params.Q() = mat_type::Zero(4,4);
    params.P() = mat_type::Zero(4,4);
    params.Q()(1,1) = 10; params.Q()(3,3) = 0;
    params.P()(0,0) = 10;
    params.R()(0,0) = 1; params.R()(1,1) = 1;
    params.x_des = []( const double & /*t*/, const state_type & /*x*/,
		       vec_type & xdes) { xdes << 1.0, 1.0, 0, 0; };
    params.mxdes_tf() << 1.0, 1.0, 0, 0;

    /* initialize SAC stepper */
    sac_step SACit( params );

    /* simulation start/stop times */
    double t0=0.0, tsim = 10;
    
    /* initial state and control */
    state_type x0(params.xlen());    x0 = { 0, 0.5, 0, 0 };
    b_control u1(params);     u1.stimes( 0, params.calc_tm() );

    /* compute final trajectory cost */
    vector<state_type> x_out, u2list, TiTappTf;
    vector<double> state_times;
    state_intp traj_intp(x_out, state_times, params.xlen());
    traj_cost J_traj( traj_intp, u2list, TiTappTf, params );
    //
    J_traj.Q() = params.Q(); J_traj.P() = params.P();
    J_traj.R() = params.R();
    //]
    
    /* Save initial x and u */
    x_out.push_back(x0); state_times.push_back(t0);
    u2list.push_back( u1.m_u_switch ); 
    TiTappTf.push_back( { u1.m_tau1, 0.0, u1.m_tau2 } );

    // [ TESTING
    // cout << "testing simHybX\n";
    // state_type u_(2), x_(4); 
    // b_control uu(params);
    // sys_dynam xdot_(uu);
    // vector<state_type> x_vec;
    // vector<double> times, events;
    // x_ = {0.592276, 0.729348, 1.12892, 1.47247};
    // u_ = { -0.700937, 16.311 };
    // uu = u_;  uu.stimes(0.52, 0.545);
    // simHybX( xdot_, x_, 0.545, 1.11, x_vec, times, events );
    // cout << "test complete\n";
    //]

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



