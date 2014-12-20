#ifndef INC_COST_HPP
#define INC_COST_HPP

namespace sac {

  //! \warning Class MUST BE MODIFIED BY USER to accomodate angle wrapping
  /*!
    General incremental trajectory tracking cost, \f$l(x)\f$, for integration
    \f$J_1 = \int_{t_0}^{t_f} l(x) \, dt + m(x(t_f))\f$.  Keeps references to 
    state interpolator so that changes state trajectory are automatically 
    accounted for.
  */
  class inc_cost {
    Eigen::Matrix< double, xlen, 1 > mx_;
    void (*p_get_DesTraj)( const double t,        // pointer to function
			   Eigen::Matrix< double, // to get desired trajectory
			   xlen, 1 > &mxdes );
  
  public:
    state_intp & m_x_intp; // store current state
    state_type m_x;
    Eigen::Matrix< double, xlen, 1 > m_mxdes;
  
    //! \todo Alex: make inputs const ref type
    /*!
      Initializes references to user maintained trajectory object and a pointer
      to the desired trajectory.
      \param[in] x_intp User maintained state interpolation object
      \param[in] xdesFnptr Pointer to the desired trajectory
    */
    inc_cost( state_intp & x_intp,
	      void (* xdesFnptr) ( const double t, 
				   Eigen::Matrix< double, xlen, 1 > &mxdes ) 
	      ): mx_(Eigen::Matrix< double,xlen,1 >::Zero(xlen,1)), 
		 p_get_DesTraj( xdesFnptr ),
		 m_x_intp( x_intp ), m_x( xlen ),
		 m_mxdes(Eigen::Matrix< double,xlen,1 >::Zero(xlen,1)) { }
  
    /*!
      Computes the value of incremental trajectory tracking cost, \f$l(x)\f$.
      \param[in] J The current cost
      \param[out] dJdt The previous incremental cost
      \param[in] t The current time
    */
    void operator() (const state_type &/*J*/, state_type &dJdt, const double t)
    {
      m_x_intp(t, m_x); // store the current state in x
      State2Mat( m_x, mx_ ); // convert state to matrix form
      //
      p_get_DesTraj( t, m_mxdes ); // Store desired trajectory point in m_mxdes
      //
      dJdt[0] = ( ( (mx_-m_mxdes).transpose() * Q 
		    * (mx_-m_mxdes) ) / 2.0 )[0];
    }

    /*!
      Computes the value of the derivative of the incremental cost, 
      \f$D_x l(x)\f$.
      \param[in] t The current time
      \param[in] mx The current state
      \param[out] dldx The derivative \f$D_x l(x)\f$.
    */
    inline void dx( const double t, const Eigen::Matrix< double, xlen, 1 > &mx,
		    Eigen::Matrix< double, 1, xlen > &dldx ) { 
      p_get_DesTraj( t, m_mxdes ); // Store desired trajectory point in m_mxdes
      //
      dldx = (mx-m_mxdes).transpose()*Q;
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
