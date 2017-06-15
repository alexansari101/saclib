#ifndef SYS_DYNAM_HPP
#define SYS_DYNAM_HPP

namespace sac {

  //[ The rhs of x' = f(x) defined as a class
  // USER SPECIFIED:
  class sys_dynam {
    double g_;
    double h_;
    b_control & u_;
    state_type u_curr_;
  
  public:
    sys_dynam( b_control & u ) : g_(9.81), h_(2.0), u_(u), 
				 u_curr_(1) { }
  
    void operator() (const state_type &x, state_type &dxdt, const double t)
    {
      u_(t, u_curr_);
      dxdt[0] = x[1];
      dxdt[1] = ( g_*sin(x[0]) + u_curr_[0]*cos(x[0]) )/h_;
    }
  };
  //]

}

#endif  // SYS_DYNAM_HPP
