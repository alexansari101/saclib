#ifndef B_CONTROL_HPP
#define B_CONTROL_HPP

namespace sac {

  /*! 
    This class stores switching control vectors.  The control switches from
    a default, nominal control to the switching control signal when \f$\tau_1 
    \leq t \leq \tau_2\f$.
  */
  class b_control {
    state_type u_default_;
  public:
    size_t m_len;
    double m_tau1, m_tau2;
    state_type m_u_switch;
  
    //! \todo Alex: Make constructor explicit if it takes single arg.
    /*!
      Constructor sets the nominal control to the zero vector.
      \param[in] len The length of the control vector.
    */
    b_control( const size_t len ) : u_default_(len), m_len(len), 
				    m_tau1(0.0), m_tau2(0.0), 
				    m_u_switch(len) { 
      for ( size_t i=0; i<m_len; i++ ) {
	u_default_[i] = 0.0;  // set default control to 0 vector
	m_u_switch[i] = 0.0;
      }
    }
  
    /*!
      Returns the control at time t.
      \param[in] t The time at which to get the current control.
      \param[out] u_curr The control vector at time t.
    */
    void operator() ( const double t, state_type& u_curr )
    {
      if ( (t >= m_tau1) && (t <= m_tau2) ) { u_curr = m_u_switch; } 
      else { u_curr =  u_default_; }
    }

    //! \todo Alex: Think about using getter/setter functions to set u_switch
    //!  when saturation is desired and directly accessing u_switch when not.
    /*!
      Sets the value of the switching control when \f$\tau_1 \leq t \leq 
      \tau_2\f$.  Also applies saturation to the vector.
      \param[in] u_switch The desired value of the switching control when 
      \f$\tau_1  \leq t \leq \tau_2\f$.
    */
    void operator= ( const state_type & u_switch ) { 
      double val;
      for ( size_t i=0; i<m_len; i++ ) { // saturate controls
	val = u_switch[i];
	if ( val > usat[i][0] ) { val = usat[i][0]; }
	else if ( val < usat[i][1] ) { val = usat[i][1]; }
	m_u_switch[i] = val;
      }
    }

    /*!
      Sets the value of the switching control when \f$\tau_1 \leq t \leq 
      \tau_2\f$ without applying saturation.
      \param[in] u_switch The desired value of the switching control when 
      \f$\tau_1  \leq t \leq \tau_2\f$.
    */
    void no_saturate( const state_type & u_switch ) { m_u_switch = u_switch; }

    /*!
      Sets switching times \f$\tau_1\f$ and \f$\tau_2\f$ where the switching
      control is applied when \f$\tau_1 \leq t \leq \tau_2\f$.
      \param[in] t1, t2 Desired values of switching times \f$\tau_1\f$ and 
      \f$\tau_2\f$.
    */
    void stimes( const double t1, const double t2 ) { m_tau1 = t1; m_tau2 = t2; }

    /*!
      Re-set the switching control to the zero vector.  Sets switching times 
      \f$\tau_1\f$ and \f$\tau_2\f$ to 0.
    */
    void clear( ) {
      for ( size_t i=0; i<m_len; i++ ) {
	m_u_switch[i] = 0.0;
      }
      m_tau1 = 0.0;    m_tau2 = m_tau1;
    }
  };

}

#endif // B_CONTROL_HPP
