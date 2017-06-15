/*
 Full state cart pendulum system under acceleration control
 */

/******************************************************************/
/* Reserve max matrix NxN dims if known at compile time for speed */
/* Must at least be the max(x_len,u_len).                         */
#define CONST_MAT 2

#include <iostream>
#include <master.hpp>               // Master include file
#include <boost/timer/timer.hpp>    // timing stuff 
#include <traj_cost.hpp>            // Add to compute final trajectory cost

using namespace sac;

void init_test( Params &params, std::vector<std::vector<double> > range );

int main(int /* argc */ , char** /* argv */ )
{
    using namespace std;

    boost::timer::auto_cpu_timer t(4, "%w seconds\n");

    /* initialize SAC parameters */
    Params params(2,1);
    params.T() = 0.5;
    params.lam() = -5;
    params.maxdt() = 0.2;
    params.ts() = 0.01;
    params.usat() = { {10, -10} };
    params.calc_tm() = 0; // params.ts();
    params.u2search() = true;
    params.Q() = mat_type::Zero(2,2);
    params.P() = mat_type::Zero(2,2);
    // params.Q()(0,0) = 1000; params.Q()(1,1) = 10;
    params.R() << 1;
    params.proj = []( state_type & x ) { AngleWrap( x[0] ); };

    /* initialize SAC stepper */
    sac_step SACit( params );

    /* simulation start/stop times */
    double t0=0.0, tsim = 9;
    
    /* initial state and control */
    state_type x0(params.xlen());    x0 = { PI, 0 };
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

    while ( t0 < tsim ) {
      
      if (t0 < params.ts()) { params.u2search() = false; }
      
      /* Perform SAC step */
      SACit( t0, x0, u1 ); // update: u1 and x0 to time t0+ts

      /* Save updated x and u */
      x_out.push_back(x0); state_times.push_back(t0);
      u2list.push_back( u1.m_u_switch ); 
      TiTappTf.push_back( { SACit.t_i, SACit.t_app, SACit.t_f } );

      params.u2search() = true;

    } /* end main while() loop */

    t.stop(); t.report();

    /* Save to file */
    SaveVec( x_out, "x.csv" );   SaveVec( u2list, "u2list.csv" );
    SaveVec( TiTappTf, "TiTappTf.csv" );

    //[ /* Compute final trajectory cost */
    J_traj.compute_cost( 0.0, state_times.back() );
    J_traj.print();
    //]    


    // /* For State-Space Control Density Plots */
    // vector<state_type> x0uVec;
    // state_type x0u(3),x0curr(2);
    // double dt=0.05;//params.ts();
    
    // for ( x0[0] = -PI/2.0 ; x0[0] < 5.0*PI/2.0 ; x0[0]=x0[0]+dt ) {
    //   for ( x0[1] = -5 ; x0[1] < 5 ; x0[1]=x0[1]+dt ) {

    // 	x0u[0] = x0[0]; x0u[1] = x0[1];
    // 	u1.clear(); x0curr = x0;
    // 	/* Perform SAC step */
    // 	SACit( t0, x0curr, u1 ); // update: u1 and x0 to time t0+ts
 
    // 	/* update initial states, time horizon, & save control */
    // 	  if ( u1.m_tau2 - u1.m_tau1 > 0.000001 ) {
    // 	    x0u[2] = (u1.m_u_switch)[0];
    // 	  } else { x0u[2] = 0.0;}

    // 	x0uVec.push_back( x0u );

    //   } /* end main while() loop */
    // }

    // t.stop(); t.report();

    // SaveVec( x0uVec, "x0_u.csv");




    // /* INITIAL CONDITION TESTS */
    // vector<vector<double> > range(2); /* { { min, max, dx }, ... } */
    // range[0].resize(3);  range[0]= { -PI+.1, PI-.1, .1 };
    // range[1].resize(3);  range[1]= { 0, .05, .1 };
    // //
    // // init_test( params, range );
    
  
    return 0;
}



void init_test( Params &params, std::vector<std::vector<double> > range ) {

  using namespace std;
  using namespace boost::numeric::odeint;

  double max_sim_time=25;

  std::ofstream out("results.txt");
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old output buf
  std::cout.rdbuf(out.rdbuf()); //redirect std::cout to results file

  cout << "Testing over range: { {"
       << range[0][0] << ", " << range[0][1] << ", " <<range[0][2] << "}, {"
       << range[1][0] << ", " << range[1][1] << ", " <<range[1][2] 
       << "} }\n\n";

  cout << "Testing params:\n\n";
  cout << "Q = \n" << params.Q() << "\n";
  cout << "P = \n" << params.P() << "\n";
  cout << "R = \n" << params.R() << "\n";
  cout << "T = " << params.T() << "\n";
  cout << "lam = " << params.lam() << "\n";
  cout << "ts = " << params.ts() << "\n";
  cout << "calc_tm = " << params.calc_tm() << "\n";
  cout << "u2search = " << params.u2search() << "\n";
  cout << "usat = {" << params.usat()[0][0] << ", "
       << params.usat()[0][1] << "}\n\n";

  cout << "max allowed simulation time for convergence = " 
       << max_sim_time << "\n";

  boost::timer::auto_cpu_timer t(4, "%w seconds\n");

  // INITIAL CONDITION TESTS
  bool success = false; bool u2search = params.u2search();
  for ( double dth=range[0][0] ; dth <= range[0][1] ; dth+=range[0][2] ) {
    for ( double ddth=range[1][0] ; ddth <= range[1][1] ; ddth+=range[1][2] ) {
      /* initialize SAC stepper */
      sac_step SACit( params );

      /* simulation start/stop times */
      double t0=0.0, tsim = max_sim_time;
    
      /* initial state and control */
      state_type x0(params.xlen());    x0 = { PI + dth, 0 + ddth };
      b_control u1(params);     u1.stimes( 0, params.calc_tm() );

      state_type xwrap(params.xlen());

      cout << "testing x0 = (";
      for (auto &i : x0)
	cout << i << ", ";
      cout << ")\t" << "dth = " << dth << "; ddth = " << ddth << "\n";

      boost::timer::auto_cpu_timer iterTmr(4, "%w seconds\n");

      while ( t0 < tsim ) {
      
	/* If not close enough for LQR stabilization */
	xwrap = x0;  AngleWrap(xwrap[0]);

	if ( std::abs(xwrap[0]) < .01 && std::abs(xwrap[1]) < .05 ) {
	  cout << "converged at t0 = " << t0 << endl;
	  success = true;
	  break;
	}
	
	if (t0 < params.ts()) { params.u2search() = false; }

	/* Perform SAC step */
	SACit( t0, x0, u1 ); // update: u1 and x0 to time t0+ts

	params.u2search() = u2search;
 
      } /* end main while() loop */

      if (success) {
	success = false;
      } else { cout << "\nfailed!!\n\n"; }

    } /* end inner for() loop */
  } /* end outer for() loop */

  cout << "Total time:\n";
  t.stop(); t.report();

  std::cout.rdbuf(coutbuf); //reset to standard output again
  
}
