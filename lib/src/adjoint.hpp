#ifndef ADJOINT_HPP
#define ADJOINT_HPP

namespace sac {

  //! \warning Class MUST BE MODIFIED BY USER to accomodate angle wrapping
  /*!
    Evaluates the rhs of \f$\dot \rho(t) = -\frac{\partial l}{\partial x}^T 
    - \frac{\partial f}{\partial x}^T \rho(t)\f$ for integration of co-state
    variable, \f$\rho(t)\f$.  Keeps references to state interpolator so that 
    changes state trajectory are automatically accounted for.
  */
  class adjoint {
    state_type x_, u1_, rho_;
    Eigen::Matrix< double, xlen, 1 > mx_, mrho_, mrhodot_;
    Eigen::Matrix< double, 1, xlen > mdldx_;
    Eigen::Matrix< double, xlen, xlen > mdfdx_;
    size_t indx_;

  public:
    state_intp & m_x_intp;
    sys_lin m_lin;
    inc_cost & m_lofx;
  
    //! \todo Alex: make inputs const ref type.
    /*!
      Initializes references to user maintained state interpolation object and a
      cost object.
      \param[in] x_intp User maintained state interpolation object
      \param[in] J The cost object required to provide incremental cost partial, 
      \f$\frac{\partial l}{\partial x}\f$
    */
    adjoint( state_intp & x_intp,
	     cost & J ) :  x_(xlen), u1_(ulen),
			   rho_(xlen), m_x_intp( x_intp ),
			   m_lofx(J.m_lofx) {  
      for ( size_t i=0; i<ulen; i++ ) { u1_[i] = 0.0; } 
    }

    /*!
      Returns the rhs of \f$\dot \rho(t) = -\frac{\partial l}{\partial x}^T 
      - \frac{\partial f}{\partial x}^T \rho(t)\f$.  The dynamics of co-state
      variable, \f$\rho(t)\f$.
      \param[in] rho The co-state variable at time t
      \param[out] rhodot The dynamics of the co-state at time t
      \param[in] t The time variable
    */
    void operator() (const state_type &rho, state_type &rhodot, const double t)
    {
      rho_ = rho;
      m_x_intp(t, x_);        // store the current state in x
      State2Mat( x_, mx_ );   // convert state to matrix form
      State2Mat( rho_, mrho_ );
      //
      m_lin.A( x_, u1_, mdfdx_ );
      //
      m_lofx.dx( t, mx_, mdldx_ );
      //
      mrhodot_ = -mdldx_.transpose() - mdfdx_.transpose()*mrho_;
      //
      for (indx_ = 0; indx_ < xlen; indx_++ ) { rhodot[indx_] = mrhodot_[indx_]; }
    }
  };

}

#endif  // ADJOINT_HPP
