#ifndef SYS_DYNAM_HPP
#define SYS_DYNAM_HPP

namespace sac {

  const size_t xlen = 4;
  const size_t ulen = 1;
  const double p1 = 0.0308, p2 = 0.0106, p3 = 0.0095,
    p4 = 0.2087, p5 = 0.0629;


  //[ The rhs of x' = f(x) defined as a class
  // USER SPECIFIED:
  class sys_dynam {
    double g_;
    b_control & u_;
    state_type u_curr_;
  
  public:
    sys_dynam( b_control & u ) : g_(9.81) , u_(u) , u_curr_(ulen) { }
  
    void operator() (const state_type &x, state_type &dxdt, const double t)
    {
      u_(t, u_curr_);
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
    }
  };
  //]

}

#endif  // SYS_DYNAM_HPP
