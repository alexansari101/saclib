#ifndef INC_COST_HPP
#define INC_COST_HPP

namespace sac {

  //! \warning Projection and derivative of projection in testing
  /*!
    General incremental trajectory tracking cost, \f$l(x)\f$, for integration
    \f$J_1 = \int_{t_0}^{t_f} l(x) \, dt + m(x(t_f))\f$.  Keeps references to 
    state interpolator so that changes state trajectory are automatically 
    accounted for.
  */
  class inc_cost {
    vec_type mx_;
    mat_type mgproj_x_;
    Params & p_;
  
  public:
    state_intp & m_x_intp; // store current state
    state_type m_x;
    vec_type m_mxdes;
  
    //! \todo Alex: make inputs const ref type
    /*!
      Initializes references to user maintained trajectory object and a pointer
      to the desired trajectory.
      \param[in] x_intp User maintained state interpolation object
      \param[in] p SAC parameters
    */
    inc_cost( state_intp & x_intp, Params & p )
      : mx_(vec_type::Zero(p.xlen(),1)), 
	mgproj_x_( mat_type::Identity(p.xlen(),p.xlen()) ),
	p_(p), m_x_intp( x_intp ), m_x( p.xlen() ),
	m_mxdes( vec_type::Zero(p.xlen(),1) ) { }
  
    /*!
      Computes the value of incremental trajectory tracking cost, \f$l(x)\f$.
      \param[in] J The current cost
      \param[out] dJdt The previous incremental cost
      \param[in] t The current time
    */
    void operator() (const state_type &/*J*/, state_type &dJdt, const double t)
    {
      m_x_intp(t, m_x); // store the current state in x
      p_.proj( m_x );
      //
      p_.x_des( t, m_x, m_mxdes ); // Get desired trajectory point
      //
      State2Mat( m_x, mx_ ); // convert state to matrix form
      //
      dJdt[0] = ( ( (mx_-m_mxdes).transpose() * p_.Q() 
      		    * (mx_-m_mxdes) ) / 2.0 )(0);

      // dJdt[0] = p_.inc_x_cost( t, m_x );
    }

    //! \todo: Alex: is it ok that this needs to be called with proj(mx)?
    /*!
      Computes the value of the gradient of the incremental cost, 
      \f$D_x l(x)^T\f$.
      \param[in] t The current time
      \param[in] mx The current state
      \param[out] dldx The gradient \f$D_x l(x)^T\f$.
    */
    inline void grad( const double t, const vec_type &mx,
		      vec_type &Glofx ) { 
      Mat2State( mx, m_x);
      p_.gproj( m_x, mgproj_x_ );
      //
      p_.x_des( t, m_x, m_mxdes ); // Get desired trajectory point
      //
      Glofx = mgproj_x_*p_.Q()*(mx-m_mxdes);

      // p_.dinc_x_cost( t, mx, dldx );
    }

    /*!
      Returns the initial time for integration, \f$t_0\f$
      \return Initial integration time \f$t_0\f$
    */
    double begin( ) { return m_x_intp.begin( ); }

    /*!
      Returns the final time for integration, \f$t_f\f$
      \return Final integration time \f$t_f\f$
    */
    double end( ) { return m_x_intp.end( ); }
  };

}

#endif  // INC_COST_HPP
