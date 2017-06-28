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


#ifndef U2_COST_HPP
#define U2_COST_HPP

namespace sac {

  /*! 
    Provides a cost function of time that can be optimized to find the best time 
    to apply the optimal control law u2*. e.g. \f$cost(t) = \sqrt{u_2^*(t)^T 
    u_2^*(t)} + \frac{dJ_1}{d \lambda^+}(t) + t^{1.6}\f$
  */
  class u2_cost {
    u2_optimal & u2Opt_;
    mode_insert_grad & dJdlam_;
    state_type u2Opt_curr_;
    double dJdlam_curr_;
    size_t i_;
    double u2_sum_;
    Params &p_;
  
  public:
    //! \todo Alex: see about converting inputs to const ref types.
    /*! 
      Initializes references to user defined objects \f$u_2^*(t)\f$ and 
      \f$\frac{dJ_1}{d \lambda^+}(t)\f$ which automatically update along with
      updates in trajectory.
      \param[in] u2Opt The object storing the values of control \f$u_2^*(t)\f$
      \param[in] dJdlam The object storing the values of mode insertion gradient 
      \f$\frac{dJ_1}{d \lambda^+}(t)\f$
      \param[in] p SAC parameters
    */
    u2_cost( u2_optimal & u2Opt, mode_insert_grad & dJdlam,
	     Params &p ) : u2Opt_( u2Opt ) , 
			   dJdlam_( dJdlam ),
			   u2Opt_curr_(p.ulen()), 
			   dJdlam_curr_(0.0), p_(p) { }
  
    /*! 
      Returns the value of a cost function at the specified time. The cost 
      function can be searched to find the best time to apply \f$u_2^*(t)\f$.
      \f[cost(t) = \sqrt{u_2^*(t)^T u_2^*(t)} + \frac{dJ_1}{d \lambda^+}(t) 
      + t^{1.6}\f]
      \param[in] t The time at which to evaluate the cost
      \return The value of the cost at the specified time
    */
    double operator() ( const double t )
    {
      u2Opt_(t, u2Opt_curr_);
      u2_sum_ = 0;    
      for ( i_ = 0; i_ < p_.ulen(); i_++ ) {u2_sum_ +=pow(u2Opt_curr_[i_], 2);}
      u2_sum_ = pow(u2_sum_, 0.5);
      //
      dJdlam_(t, dJdlam_curr_);
      //
      return u2_sum_ + dJdlam_curr_ + pow(t, 1.6);
    }
  };

}

#endif  // U2_COST_HPP
