/*
 Pendubot system under torque control
 */

/******************************************************************/
/* Reserve max matrix NxN dims if known at compile time for speed */
/* Must at least be the max(x_len,u_len).                         */
#define CONST_MAT 4

#include <iostream>
#include <boost/timer/timer.hpp>           // timing stuff 
#include <master.hpp>                      // Master include file
#include <traj_cost.hpp>                   // Add to compute final trajectory cost

using namespace sac;


//[ The rhs of x' = f(x,-kx) defined as a class
// USER SPECIFIED:
class fb_sys_dynam {
  double g_;
  state_type u_curr_;
  
public:
  double ufb;

  fb_sys_dynam( ) : g_(9.81) , u_curr_(ulen) { }
  
  void operator() (const state_type &x, state_type &dxdt, const double /*t*/)
  {
    u_curr_[0] = ufb;
    //
    dxdt[0] = x[1];
    dxdt[1] = ( g_*(2*p2*p4 - p3*p5)*sin(x[0]) + 
		g_*p3*p5*sin(x[0]-2*x[2]) + 2*p2*u_curr_[0] - 
		2*p3*sin(x[0]-x[2])*(p3*cos(x[0]-x[2])*pow(x[1], 2) + 
				     p2*pow(x[3], 2))
		)/( 2*p1*p2 - 2*pow(p3, 2)*pow(cos(x[0]-x[2]), 2) );
    dxdt[2] = x[3];
    dxdt[3] = (-g_*p3*p4*cos(x[0]-x[2])*sin(x[0]) + g_*p1*p5*sin(x[2]) -
	       p3*cos(x[0]-x[2])*u_curr_[0] + p3*sin(x[0]-x[2])*
	       (p1*pow(x[1], 2) + p3*cos(x[0]-x[2])*pow(x[3], 2))
	       )/( p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2) );
    //
  }
};
//]




int main(int /* argc */ , char** /* argv */ )
{
    using namespace std;
    using namespace boost::numeric::odeint;

    boost::timer::auto_cpu_timer t(4, "%w seconds\n");

    /* initialize SAC parameters */
    Params params(4,1);
    params.T() = 0.5;
    params.lam() = -15;
    params.maxdt() = 0.2;
    params.ts() = 0.005;
    params.usat() = { {7, -7} };
    params.calc_tm() = params.ts();
    params.u2search() = true;
    params.Q() = mat_type::Zero(4,4);
    params.P() = mat_type::Zero(4,4);
    params.Q()(0,0) = 100; params.Q()(1,1) = 0.0001; 
    params.Q()(2,2) = 200; params.Q()(3,3) = 0.0001;
    params.R() << 0.1;
    params.proj = []( state_type & x ) { AngleWrap( x[0] ); 
					 AngleWrap( x[2] ); };

    /* initialize SAC stepper */
    sac_step SACit( params );

    /* simulation start/stop times */
    double t0=0.0, tsim = 20;
    
    /* initial state and control */
    state_type x0(params.xlen());    x0 = { PI, 0, PI, 0 };
    b_control u1(params);     u1.stimes( 0, params.calc_tm() );

    /* LQR stabilization */
    fb_sys_dynam xdotfb;  state_type xwrap(xlen);
    runge_kutta4< state_type > stepper;
    bool stabilizing = false; size_t count = 0;

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
      
      /* If not close enough for LQR stabilization */
      xwrap = x0;  AngleWrap(xwrap[0]); AngleWrap(xwrap[2]);
      if ( !stabilizing &&
	   ( std::abs(xwrap[0]) > .05 /*|| std::abs(xwrap[1]) > .25*/
	     || std::abs(xwrap[2]) > .05 /*|| std::abs(xwrap[3]) > .25*/ ) ) {

	/* Perform SAC step */
	SACit( t0, x0, u1 ); // update: u1 and x0 to time t0+ts

	/* Save updated x and u */
	x_out.push_back(x0); state_times.push_back(t0);
	u2list.push_back( u1.m_u_switch ); 
	TiTappTf.push_back( { SACit.t_i, SACit.t_app, SACit.t_f } );

      }
      else {   /* else do LQR stabilization */

	if ( count < 1 ) { 
	  cout << "stabilizing at t = " << t0 << "\n"; 
	  (TiTappTf.back()).back() = t0;
	  count++; 
	  stabilizing = true; }
	
	/* Use the same feedback [t0, t0+ts] */
	xdotfb.ufb= -(-0.231168*xwrap[0]-1.73764*xwrap[1]
		      -28.9858*xwrap[2]-3.85553*xwrap[3]);
	
	/* Simulate LQR feedback stabilized trajectory */
	integrate( xdotfb , x0 , t0 , t0+params.ts(), params.ts() );

	t0=t0+params.ts();

	/* Save updated x and u */
	x_out.push_back(x0); state_times.push_back(t0);
	u2list.push_back( { xdotfb.ufb } );
	TiTappTf.push_back( { t0, t0, t0+params.ts() } );
      }
 
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
