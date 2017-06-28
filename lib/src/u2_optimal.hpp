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


#ifndef U2_OPTIMAL_HPP
#define U2_OPTIMAL_HPP

namespace sac {

  /*! 
    Stores the optimal switching control,\f$u_2^*(t)\f$.  Keeps references to 
    user maintained state and co-state trajectory interpolation objects so that 
    \f$u_2^*(t)\f$ automatically updates along with changes in trajectory.
    \f[u_2^*(t) = \; (\Lambda + R)^{-1} \, [\Lambda \, u_1(t) + h(x(t))^T \rho(t)
    \, \alpha_d]\f]
    where \f$\Lambda = h(x(t))^T \rho(t) \rho(t)^T h(x(t))\f$ and \f$h(x(t)) = 
    \frac{\partial f_1}{\partial u_1}\f$.  The dynamics, \f$f_1\f$, should be in
    control affine form.
  */
  class u2_optimal {
    state_intp & rx_intp_;
    state_intp & rrho_intp_;
    state_type x_curr_;
    state_type rho_curr_;
    double & alpha_;
    state_type u1_;
    vec_type  mrslt_, mu1_;
    vec_type mx_curr_, mrho_curr_;
    mat_type B_;
    sys_lin lin_;
    size_t i_;
    Params & p_;
  
  public:
    //! \todo Alex: make inputs const ref type
    /*!
      Initializes references to user maintained trajectory objects and the 
      desired rate of change of a trajectory tracking cost functional.
      \param[in] x_intp User maintained state interpolation object
      \param[in] rho_intp User maintained co-state interpolation object
      \param[in] alpha User specified desired change in cost functional
      relative to the duration of activiation of \f$u_2^*(t)\f$. 
      i.e. \f$\frac{\Delta J_1}{\Delta t}\f$.
      \param[in] p SAC parameters
    */
    u2_optimal( state_intp & x_intp, state_intp & rho_intp, 
		double & alpha, Params & p ) : rx_intp_( x_intp ) ,
				   rrho_intp_( rho_intp ) ,
				   x_curr_(p.xlen()) , 
				   rho_curr_(p.xlen()) , 
				   alpha_( alpha ) , u1_(p.ulen()),
				   mrslt_(p.ulen(),1), mu1_(p.ulen(),1),
				   mx_curr_(p.xlen(),1), 
				   mrho_curr_(p.xlen(),1),
				   B_(p.xlen(),p.ulen()), p_(p) { 
      for ( i_=0; i_<p.ulen(); i_++ ) {
	u1_[i_] = 0.0;
	mu1_(i_,0) = 0.0;
	mrslt_(i_,0) = 0.0;
      }
    }
  
    //! \todo Alex: incorporate u1
    /*!
      Computes the value of the optimal switching control,\f$u_2^*(t)\f$.
      \f[u_2^*(t) = \; (\Lambda + R)^{-1} \, [\Lambda \, u_1(t) + h(x(t))^T 
      \rho(t) \, \alpha_d]\f]
      \param[in] t The time at which to compute the mode insertion gradient
      \param[out] u2Opt_curr The value of the optimal switching control
    */
    void operator() ( const double t, state_type & u2Opt_curr ) {
      rx_intp_(t, x_curr_);      State2Mat( x_curr_, mx_curr_ );
      rrho_intp_(t, rho_curr_);  State2Mat( rho_curr_, mrho_curr_ );
      lin_.B( x_curr_, u1_, B_ );
      mrslt_ = B_.transpose()*mrho_curr_;
      mrslt_ = ( mrslt_*mrslt_.transpose() + p_.R() ).inverse()
	*(/* (...)*mu1_ +*/ mrslt_*alpha_);    
      //
      for ( i_=0; i_<p_.ulen(); i_++ ) { u2Opt_curr[i_] = mrslt_(i_,0); }
    }

  };

}

#endif  // U2_OPTIMAL_HPP
