/*
 2D Compliant (spring-leg) locomotion

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
 Examples use k = (14 kN/m and) 20 kN/m
 * massless legs
 * point mass at COM
 * walking is around 1.2 m/s and running is around 4.0 m/s and starts > 3 m/s

 Phases: 
 Left single support -> double support -> right single support

 Switching conditions:
 Touchdown -> z_TD = L0 sin(alpha_0)
 Lift-off -> L = L0

 State dimension is 7:
 X = [ x,    x',   z,    z',   t_lx, t_rx, alpha_0, loc ];
 X = [ x[0], x[1], x[2], x[3], x[4], x[5], x[6],    x[7]];
 
 Control dimension is 1:
 U = [ dalpha_0 ];

 Left leg single support dynamics:
 fl = [ x',
        Cl*(x-t_lx)/m,
	z',
	(Cl*(z-zGrnd(t_lx)) - mg)/m,
	0,
	0,
	u[0], 
	0 ];

 Double support dynamics:
 fd = [ x',
        ( Cl*(x-t_lx) - Cr*(t_rx-x) )/m,
	z',
	( Cl*(z-zGrnd(t_lx)) + Cr*(z-zGrnd(t_rx)) - mg )/m,
	0,
	0,
	u[0],
	0 ];

 Right leg single support dynamics:
 fr = [ x',
        -Cr*(t_rx-x)/m,
	z',
	(Cr*(z-zGrnd(t_rx)) - mg)/m,
	0,
	0,
	u[0],
	0 ];

 Ll = sqrt((x-t_lx)^2 + (z-zGrnd(t_lx))^2)
 Lr = sqrt((t_rx-x)^2 + (z-zGrnd(t_rx))^2)
 Cl = k*( L0/Ll - 1)
 Cr = k*( L0/Lr - 1)
 d = FP_{i+1,x} - FP_{i,x} = t_rx - t_lx
 * FP is the foot point of a stance spring - t_r/lx is right/left toe x-coord
 * t_lx and t_rx reset impulsively at phase transitions
 */

/******************************************************************/
/* Reserve max matrix NxN dims if known at compile time for speed */
/* Must at least be the max(x_len,u_len).                         */
#define CONST_MAT 10

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
    Params params(8,2);
    params.T() = 1;
    params.lam() = -10;
    params.maxdt() = 0.2;
    params.ts() = 0.01;
    params.usat() = { {.1, -.1}, {1, -1} };
    params.calc_tm() = params.ts();
    params.u2search() = false;
    params.Q() = mat_type::Zero(params.xlen(),params.xlen());
    params.P() = mat_type::Zero(params.xlen(),params.xlen());
    params.R() = mat_type::Identity(params.ulen(),params.ulen());
    params.Q()(1,1) = 100; // params.Q()(3,3) = 100;// params.Q()(2,2) = 100;
    params.x_des = []( const double & /*t*/, const state_type & /*x*/,
    		       vec_type & xdes) { xdes << 0, 1.2, 1, 0, 
    					  0, 0, 69.0*PI/180.0, 0; };

    /* initialize SAC stepper */
    sac_step SACit( params );

    /* simulation start/stop times */
    double t0=0.0, tsim = 10;
    
    /* initial state and control */
    state_type x0(params.xlen());    x0 = { 0, 0 /* .883 for Es = 816J */, 
    					    1-.056, 0, 0, 0, 69*PI/180.0, 0 };
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

    //[ Terminal State
    cout << "x0: ";
    for (auto &i : x0)
      cout << i << ", ";
    cout << "\n";
    /**/
  
    return 0;
}



