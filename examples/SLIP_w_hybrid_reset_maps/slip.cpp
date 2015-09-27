/*
 3D SLIP Hopper

@article{DblSpringWalkRun2006,
  title={Compliant leg behaviour explains basic dynamics of walking and running},
  author={Geyer, Hartmut and Seyfarth, Andre and Blickhan, Reinhard},
  journal={Proceedings of the Royal Society of London B: Biological Sciences},
  volume={273},
  number={1603},
  pages={2861--2867},
  year={2006},
  publisher={The Royal Society}
}

% First suggests the double spring model for compliant walking (it is also used for running) and provides evidence that walking is like running in that it depends on bouncing gaits and that stiff pendulum walking models are not sufficient.  Walking and running are separated by a gap of 1.5 m/s in locomotion speed or, equivalently, a gap in system energy.  Stiff pendulum models also lack double support phase and use parameters to determine how much stride energy is recovered during walking.  The compliant walking model shows double the double support phase is critical and it is unnecessary to reduce stride energy.

 Params (from Geyer paper above):
 m = 80 kg; L0 = 1 m; g = 9.81 m/s^2
 Walking examples use k = (14 kN/m and) 20 kN/m
 * massless legs
 * point mass at COM
 * walking is around 1.2 m/s and running is around 4.0 m/s and starts > 3 m/s

 State dimension is 8:
 X =  [ x,    x',   y,    y',   z,    z',   x_toe,   y_toe,  loc ];
 X =  [ x[0], x[1], x[2], x[3], x[4], x[5], x[6],    x[7],   x[8]];
 
 Control dimension is 3 with saturation:
 U = [ add_toe_x_vel,   add_toe_y_vel,   thrust_into_ground ];

 Stance dynamics:
 fs = [ x',
        ( k/m ( L0-L(x,y,z) )  + u3 )*(x-xt)/L,
	y',
	( k/m ( L0 - L(x,y,z) ) + u3 )*(y-yt)/L,
	z',
	( k/m ( L0 - L(x,y,z) ) + u3 )*(z - zGrndToe)/L - g,
	0,
	0,
	0 ];

 Flight dynamics:
 ff = [ x', 0, y', 0, z', -g, x' + u1, y' + u2, 0 ];

 Dynamics switch when:
 z = L0*(z - zGrndToe) / Sqrt[ (x - xt)^2 + (y - yt)^2 + (z - zGrndToe)^2 ]

 Change Log:
 + added a line to SimInitXRho to filter event that get caught after tf
     I catching events after tf and it was integration causing issues
     when I tried to simulated from [t_event, tf] (it reversed time).
 + SimHybX: still having issues with high frequency switing.  I do some
   filtering before I throw to make sure that the system velocity is away
   from the guard and that z > 0.
 + SimHybX: I'm overflowing doubles (e.g., 1e-324**2). I need to check range.
 */

/******************************************************************/
/* Reserve max matrix NxN dims if known at compile time for speed */
/* Must at least be the max(x_len,u_len).                         */
#define CONST_MAT 9

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
    Params params(9,3);
    params.T() = 0.6;
    params.lam() = -10;
    params.maxdt() = 0.2;
    params.ts() = 1.0/100.0;
    params.usat() = { {5, -5}, {5, -5}, {30, -30} };
    params.calc_tm() = params.ts();
    params.u2search() = false;
    params.Q() = mat_type::Zero(params.xlen(),params.xlen());
    params.P() = mat_type::Zero(params.xlen(),params.xlen());
    params.R() = mat_type::Identity(params.ulen(),params.ulen());
    params.Q()(1,1) = 70; params.Q()(3,3) = 70; params.Q()(4,4) = 50;
    params.mxdes_tf() << 0, .7, 0, .7, 1.25, 0, 0, 0, 0;

    /* initialize SAC stepper */
    sac_step SACit( params );

    /* simulation start/stop times */
    double t0=0.0, tsim = 10;
    
    /* initial state and control */
    state_type x0(params.xlen());    x0 = { 0, 0, 0, 0, 1.5, 0, 0, 0, 3 };
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
    
    // /* Save initial x and u */
    // x_out.push_back(x0); state_times.push_back(t0);
    // u2list.push_back( u1.m_u_switch ); 
    // TiTappTf.push_back( { u1.m_tau1, 0.0, u1.m_tau2 } );


    // Randomly Sample Initial conditions
    state_type xorig(params.xlen()); xorig = { -.1, 0, 0, 0, 1.5, 0, 0, 0, 3 };
    
    // Initialize RNG with a seed
    srand( time(NULL) );
    // Sample points uniformly in an interval centered on the components 
    // of xorig
    double interval = 0.4;
    for (size_t its=0; its < 10; its++) {
      //[ Randomize Initial state around xorig
      x0 = xorig;
      t0=0.0;
      u1.stimes( 0, params.calc_tm() );
      cout << "\n\n x0: (";
      for (size_t i=0; i<x0.size(); i++) {
	if ( i==0 || i==2 )
	  x0[i] += ((double) rand( ) / (RAND_MAX))*interval  - interval/2;

	cout << x0[i] << ", ";
      }
      cout << ")\n";
      double rv1 = ((double) rand( ) / (RAND_MAX))*interval  - interval/2;
      double rv2 = ((double) rand( ) / (RAND_MAX))*interval  - interval/2;
      cout << "vx = " << 0.5 + rv1 << endl;
      cout << "vy = " << 0.7 + rv2 << endl;
      /**/    

      while ( t0 < tsim ) {
   
	params.x_des = [x0,rv1,rv2]( const double & t, const state_type &x,
			     vec_type & xdes) { xdes << 0, .5+rv1, 0, .7+rv2, 
						// 1.4,
						// .5*x0[0]+1.5/*1.4*/, 
						zGrndToe(x0[6],x0[7])+1.75,
						0, 0, 0, 0; };
   
	if (x0[4] < zGrndToe(x0[0],x0[2])) {
	  cout << "\nFell at t = " << t0 << "\n\n"; break; }

	/* Perform SAC step */
	SACit( t0, x0, u1 ); // update: u1 and x0 to time t0+ts

	// /* Save updated x and u */
	// x_out.push_back(x0); state_times.push_back(t0);
	// u2list.push_back( u1.m_u_switch ); 
	// TiTappTf.push_back( { SACit.t_i, SACit.t_app, SACit.t_f } );

      } /* end main while() loop */

      t.stop(); t.report();

      // /* Save to file */
      // SaveVec( x_out, "x.csv" );   SaveVec( u2list, "u2list.csv" );
      // SaveVec( TiTappTf, "TiTappTf.csv" );

      // //[ /* Compute final trajectory cost */
      // J_traj.compute_cost( 0.0, state_times.back() );
      // J_traj.print();
      // //]

      //[ Terminal State
      cout << "Terminal x: (";
      for (auto &i : x0)
	cout << i << ", ";
      cout << ")\n";
      /**/    

    } /* end for its */

    return 0;
}

