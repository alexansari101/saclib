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


#ifndef MODE_INSERT_GRAD_HPP
#define MODE_INSERT_GRAD_HPP

namespace sac {

  /*!
    Stores the values of the mode insertion gradient, \f$\frac{dJ1}{d\lambda^+}(t)
    = \rho(t)^T[ f(u_2^*(t))-f(u_1(t)) ]\f$.  Keeps references to user maintained
    state and co-state trajectory interpolation objects and \f$u_2^*(t)\f$ object
    so that the mode insertion graident automatically updates along with changes 
    in trajectory and \f$u_2^*(t)\f$.
  */
  class mode_insert_grad {
    state_intp & rx_intp_;
    state_intp & rrho_intp_;
    u2_optimal & u2Opt_;
    state_type x_curr_, rho_curr_, u2_curr_;
    b_control u1_, u2_;
    sys_dynam xdot1_, xdot2_;
    state_type dxdt1_, dxdt2_;
    vec_type mrho_curr_, mxdot1_ , mxdot2_;
  
  public:
    //! \todo Alex: make inputs const ref type
    /*!
      Initializes references to user maintained trajectory and control objects.
      \param[in] x_intp User maintained state interpolation object
      \param[in] rho_intp User maintained co-state interpolation object
      \param[in] u2Opt User maintained \f$u_2^*(t)\f$ object
      \param[in] p SAC parameters
    */
    mode_insert_grad( state_intp & x_intp, state_intp & rho_intp,
		      u2_optimal & u2Opt, Params & p ) 
      : rx_intp_( x_intp ) ,
	rrho_intp_( rho_intp ) , 
	u2Opt_(u2Opt) , x_curr_(p.xlen()) , 
	rho_curr_(p.xlen()) , u2_curr_(p.ulen()) ,
	u1_(p) , u2_(p) ,
	xdot1_(u1_) , xdot2_(u2_) , 
	dxdt1_(p.xlen()) ,  dxdt2_(p.xlen()),
	mrho_curr_(p.xlen(),1), 
	mxdot1_(p.xlen(),1) ,
	mxdot2_(p.xlen(),1) { }
  
    /*!
      Computes the value of the mode insertion gradient, \f$\frac{dJ1}
      {d\lambda^+}(t)\f$.
      \f[\frac{dJ1}{d\lambda^+}(t) = \rho(t)^T[ f(u_2^*(t))-f(u_1(t)) ]\f]
      \param[in] t The time at which to compute the mode insertion gradient
      \param[out] dJdlam_curr The value of the mode insertion gradient
    */
    void operator() ( const double t, double & dJdlam_curr )
    {
      rx_intp_(t, x_curr_);
      rrho_intp_(t, rho_curr_);
      u2Opt_( t, u2_curr_ );
      //
      u2_.no_saturate( u2_curr_ ); 
      u2_.stimes( rx_intp_.begin(), rx_intp_.end() );
      State2Mat( rho_curr_, mrho_curr_ );
      //
      xdot1_( x_curr_, dxdt1_, t );
      State2Mat( dxdt1_, mxdot1_ );
      //
      xdot2_( x_curr_, dxdt2_, t );
      State2Mat( dxdt2_, mxdot2_ );
      //
      dJdlam_curr = ( mrho_curr_.transpose()*( mxdot2_ - mxdot1_ ) )(0);
    }
  };

}

#endif  // MODE_INSERT_GRAD_HPP
