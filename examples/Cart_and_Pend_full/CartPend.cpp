/*
 Full state cart pendulum system under acceleration control
 */

#include "CartPend.hpp"

#include <boost/timer/timer.hpp>    // timing stuff 
#include <traj_cost.hpp>            // Add to compute final trajectory cost

using namespace sac;

void output( std::vector<state_type>& states, std::vector<double>& times );

int main(int /* argc */ , char** /* argv */ )
{
    using namespace std;
    using namespace sac::init;

    boost::timer::auto_cpu_timer t(4, "%w seconds\n");

    /* initialization */
    // Params params(4,1);
    params.T() = 0.28;
    params.lam() = -10;
    params.maxdt() = 0.2;
    params.ts() = 0.001;
    params.usat() = { {25, -25} };
    params.calc_tm() = ts;
    params.u2search() = false;
    params.Q() = Eigen::Matrix<double, 4,4>::Zero(4,4);
    params.P() = Eigen::Matrix<double, 4,4>::Zero(4,4);
    params.P()(0,0) = 500;
    params.R() << 0.3;

    /* simulation start/stop times */
    double t0=0.0, tsim = 10;
    
    /* initial state and control */
    state_type x0(xlen), u_default(ulen);
    x0 = { PI, 0, 0, 0 };  t_curr = { 0.0, 0.0, params.calc_tm() };
    u_switch = { 0 };

    /* compute final trajectory cost */
    vector<double> state_times;
    state_intp traj_intp(x_out, state_times, xlen);
    traj_cost J_traj( traj_intp, state_proj, u2list, TiTappTf, 
		      get_DesTraj, SACit.m_mxdes_tf );
    J_traj.Q() = 0*J_traj.Q();
    J_traj.Q()(0,0) = 1000;  J_traj.Q()(1,1) = 10;
    J_traj.P() = 0*J_traj.P();
    J_traj.R() << 0.3;
    //]
    
    /* Save initial x and u */
    x_out.push_back(x0); state_times.push_back(t0);
    u2list.push_back(u_default); TiTappTf.push_back( t_curr );

    while ( t0 < tsim ) {
      
      /* Perform SAC step */
      std::vector<double> rvec = sac_stepper(x0, u_switch, t0);
      
      /* Save initial x and u */
      x_out.push_back(x0); state_times.push_back(t0);
      u2list.push_back( u_switch ); TiTappTf.push_back( t_curr );

      t0 = t0+ts;

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
   output( rho_vec, rho_times );
   cout << "rho steps: " << steps << endl; */
void output( std::vector<state_type>& states, std::vector<double>& times ) {
  using namespace std;

  for( size_t i=0; i<times.size(); i++ )
  {
      cout << times[i] << '\t' << states[i][0] << '\t' << states[i][1] << '\n';
  }
}
