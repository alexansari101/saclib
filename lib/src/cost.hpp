#ifndef COST_HPP
#define COST_HPP

namespace sac {

  /*! 
    Class keeps track of the trajectory tracking cost.  It stores the current
    state through a reference to a state interpolator object and thus only
    requires calling an update method to re-compute trajectory cost \f$J_1 = 
    \int_{t_0}^{t_f} l(x(t))dt + m(x(t_f))\f$ after state / control updates.
  */
  class cost {
    state_type J1_, x_tf_;  // init w/ terminal cost then holds current cost
    double t0_, tf_;
    size_t J1steps_;
    Eigen::MatrixXd mx_tf_;
    Eigen::MatrixXd & mxdes_tf_;
  
  public:
    inc_cost m_lofx;          // incremental trajectory cost
    state_intp & m_x_intp;

    //! \todo Alex: Make constructor reference and pointer inputs const.
    /*!
      Constructs a cost object from a state interpolation object and desired
      state trajectory.
      \param[in] x_intp state interpolation object
      \param[in] xdesFnptr Pointer to a function that evaluates \f$x_{des}(t)\f$.
    */
    cost( state_intp & x_intp,
	  void (*xdesFnptr) ( const double t, const state_type &x,
			      Eigen::MatrixXd &m_mxdes ),
	  Eigen::MatrixXd & mxdes_tf 
	  ) : J1_(1), x_tf_( xlen ), t0_( 0.0 ), tf_( 0.0 ), 
	      J1steps_(0), mx_tf_(xlen,1), mxdes_tf_( mxdes_tf ),
	      m_lofx( x_intp, xdesFnptr ), m_x_intp( x_intp ) { }

    /*! 
      Get the cost at the terminal time, \f$m(x(t_f))\f$.
      \return \f$(x(t_f)-x_{des}(t_f))^T\;P_1\;(x(t_f)-x_{des}(t_f))\f$, 
      a quadratic form for \f$m(x(t_f))\f$.
    */
    inline double get_term_cost( ) {
      t0_ = m_lofx.begin();
      tf_ = m_lofx.end();
      m_x_intp( tf_, x_tf_ );
      State2Mat( x_tf_, mx_tf_ );
      return ((mx_tf_-mxdes_tf_).transpose()*P*(mx_tf_-mxdes_tf_))(0);
    }
  
    /*! 
      Get the derivative of the terminal time, \f$D_x m(x(t_f))\f$.
      \return \f$(x(t_f)-x_{des}(t_f))^T\;P_1\f$, which is \f$D_x m(x(t_f))\f$
      assuming a quadratic form for the terimal cost.
    */
    inline Eigen::MatrixXd get_dmdx( ) { 
      t0_ = m_lofx.begin();
      tf_ = m_lofx.end();
      m_x_intp( tf_, x_tf_ );
      State2Mat( x_tf_, mx_tf_ );
      return (mx_tf_-mxdes_tf_).transpose()*P;
    }
  
    /*!
      The function computes the integral of \f$l(x)\f$ and appends it to a 
      provided terminal cost to return cost \f$J_1 = \int_{t_0}^{t_f} l(x(t)) 
      dt + m(x(t_f))\f$.
    */
    size_t compute_cost( state_type& term_cost );
  
    /*!
      Implicitly converts cost object to a double.
      \return cost \f$J_1 = \int_{t_0}^{t_f} l(x(t))dt + m(x(t_f))\f$.
    */
    operator double() { return J1_[0]; }

    /*!
      \return Returns # of integration steps in computing cost \f$J_1 = 
      \int_{t_0}^{t_f} l(x(t))dt + m(x(t_f))\f$.
    */
    size_t steps( ) { return J1steps_; }

    /*!
      Re-computes the cost stored in the cost object.  This should be called 
      to update the cost object after state / controls have been modified.
    */
    void update( ) {
      t0_ = m_lofx.begin();
      tf_ = m_lofx.end();
      J1_[0] = get_term_cost();  // terminal cost to be added on
      J1steps_ = compute_cost( J1_ );
    }
  };


  /*!
    \param[in,out] term_cost the input terminal cost that gets updated with 
    the total cost after integration.
    \return The number of integration steps required.
  */
  size_t cost::compute_cost( state_type& term_cost ) {
    using namespace std;
    using namespace boost::numeric::odeint;
    typedef runge_kutta_dopri5< state_type > stepper_type;
    
    double eps = 1E-7;
    size_t J1_steps = integrate_adaptive( make_controlled( 1E-5 , 1E-5 , 
							   stepper_type( ) ) , 
					  m_lofx , term_cost , t0_ , tf_-eps , 0.01 );

    return J1_steps;
  }
  //]

}

#endif  // COST_HPP
