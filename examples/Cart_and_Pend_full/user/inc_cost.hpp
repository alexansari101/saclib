#ifndef INC_COST_HPP
#define INC_COST_HPP

namespace sac {

  /*********************************************/
  //[ Optimized incremental trajectory tracking cost, l(x), for integration
  // USER SPECIFIED:
  class inc_cost {
    state_intp & x_intp_; // store current state
    state_type x_;
    vec_type mxdes_;
    void (*p_get_DesTraj_)( const double t, const state_type &x,
			    vec_type &m_mxdes ); //  desired trajectory
    Params & p_;
  
  public:
    inc_cost( state_intp & x_intp, 
	      void (* xdesFnptr) ( const double t, const state_type &x,
				   vec_type &m_mxdes ),
	      Params & p
	      ) : x_intp_( x_intp ), x_( p.xlen() ),		
		  mxdes_(vec_type::Zero(p.xlen(),1) ),
		  p_get_DesTraj_( xdesFnptr ), p_(p) { }
  
    void operator() (const state_type &/*J*/, state_type &dJdt, const double t)
    {
      x_intp_(t, x_); // store the current state in x
      AngleWrap( x_[0] ); // Only for angle wrapping
      //
      p_get_DesTraj_( t, x_, mxdes_ ); // Get desired trajectory point
      //
      dJdt[0] = ( ( p_.Q()(0,0)*pow(x_[0]-mxdes_(0) , 2) 
      		    + p_.Q()(1,1)*pow(x_[1]-mxdes_(1) , 2) 
      		    + p_.Q()(2,2)*pow(x_[2]-mxdes_(2) , 2) 
      		    + p_.Q()(3,3)*pow(x_[3]-mxdes_(3) , 2)  ) / 2.0 );
    }
  
    inline void dx( const double t, const vec_type &mx,
		    mat_type &dldx ) { 
      Mat2State(mx, x_); // store the current state in x
      AngleWrap( x_[0] ); // Only for angle wrapping
      //
      p_get_DesTraj_( t, x_, mxdes_ ); // Get desired trajectory point
      //
      dldx = (mx-mxdes_).transpose()*p_.Q();
    }

    double begin( ) { return x_intp_.begin( ); }

    double end( ) { return x_intp_.end( ); }
  };
  //]

}

#endif  // INC_COST_HPP
