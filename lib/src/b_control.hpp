/**
 * A header-only library implementation of sequential action control (SAC)
 *  
 * Copyright (C) 2016 Alex R. Ansari and Todd D. Murphey
 *
 * This file is part of SAClib.
 *
 * SAClib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAClib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAClib.  If not, see <http://www.gnu.org/licenses/>.
 */


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
    Params * p_;
  public:
    double m_tau1, m_tau2;
    state_type m_u_switch;
  
    //! \todo Alex: Make constructor explicit if it takes single arg.
    //! \todo Alex: Make constructor reference and pointer inputs const.
    /*!
      Constructor sets the nominal control to the zero vector.
      \param[in] p SAC parameters
    */
    b_control( Params & p ) : u_default_(p.ulen()), p_(&p), 
				    m_tau1(0.0), m_tau2(0.0), 
				    m_u_switch(p.ulen()) { 
      for ( size_t i=0; i<p.ulen(); i++ ) {
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

    /*!
      Sets the value of the switching control when \f$\tau_1 \leq t \leq 
      \tau_2\f$.  Also applies saturation to the vector.
      \param[in] u_switch The desired value of the switching control when 
      \f$\tau_1  \leq t \leq \tau_2\f$.
    */
    void operator= ( const state_type & u_switch ) { 
      double val;
      for ( size_t i=0; i<p_->ulen(); i++ ) { // saturate controls
	val = u_switch[i];
	if ( val > p_->usat()[i][0] ) { val = p_->usat()[i][0]; }
	else if ( val < p_->usat()[i][1] ) { val = p_->usat()[i][1]; }
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
      for ( size_t i=0; i<p_->ulen(); i++ ) {
	m_u_switch[i] = 0.0;
      }
      m_tau1 = 0.0;    m_tau2 = m_tau1;
    }
  };

}

#endif // B_CONTROL_HPP
