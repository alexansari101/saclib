#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <functional>               // To use std::function

namespace sac {

  /*!
    User interface for specifying parameters required by SAC 
  */
  class Params {
    /* critical parameters */
    const size_t xlen_;   //! state dimension
    const size_t ulen_;   //! control dimension
    double T_;            //! prediction horizon
    double lam_;          //! descent rate (should be < 0)
    double maxdt_;        //! initial (max) duration for control application
    double ts_;           //! sampling time (duration between control calcs)
    std::vector< std::vector<double> > usat_;  //! control saturation ranges
    double calc_tm_;      //! time assumed for SAC calcs
    bool u2search_;       //! search for best time to apply control?

    /* cost parameters */
    mat_type mQ_;  //! incremental quadratic state cost weights 
    mat_type mP1_; //! terminal quadratic cost weights
    mat_type mR_;  //! incremental quadratic control cost weights

    /* misc parameters */
    size_t backtrack_its_; //! max # of sac backtracking iterations
    double eps_cost_;      //! cost integration error tolerance

    /* traj tracking */
    vec_type mxdes_tf_;    //! the desired trajectory at time t0 + T

  public:
    /*!
      Stores a user defined function that returns the desired trajectory
      \param[in] t time to evaluate desired trajectory
      \param[in] x state at time x(t) for evaluation of desired trajectory
      \param[out] mx_des an Eigen vec_type storing the desired trajectory
    */
    std::function<void(const double &, const state_type &, vec_type &)> x_des;

    /*!
      Projection to apply to the state matrix
      \param[out] x a state_type storing the current state
    */
    std::function<void(state_type &)> proj;

    /*!
      Derivative of the state projection
      \param[out] x a state_type storing the derivative of the projection
      at the current state
    */
    std::function<void(state_type &)> dproj;

    /*!
      Constructor for Params class to hold SAC parameters.
      \param[in] xlen state vector dimension
      \param[in] ulen control vector dimension
    */
    Params( const size_t xlen, 
	    const size_t ulen ) : xlen_(xlen), ulen_(ulen), T_(.5), 
				  lam_(-5), maxdt_(0.2), ts_(0.01), 
				  usat_(ulen, std::vector<double>(2)), 
				  calc_tm_(ts_), u2search_(false),
				  mQ_(xlen,xlen), mP1_(xlen,xlen),
				  mR_(ulen,ulen), backtrack_its_(4),
				  eps_cost_(1E-5), 
				  mxdes_tf_( vec_type::Zero(xlen,1) )
    { 
      for ( size_t i=0; i<ulen_; i++ ) {
      	usat_[i] = {100, -100}; // default saturation limits
      }
      
      x_des = [this](const double & /*t*/,const state_type & /*x*/,
		     vec_type & xdes) { xdes.Zero(this->xlen_,1); };
      proj = [](state_type & /*x*/) { };
      dproj = [](state_type & /*x*/) { };
    }

    /*!
      \return An unsigned number representing the state dimension
    */
    size_t xlen() const { return xlen_; }

    /*!
      \return An unsigned number representing the control dimension
    */
    size_t ulen() const { return ulen_; }

    /*!
      get using:      double T = param.T();
      get ref using:  double & T = param.T();
      set using:      param.T() = 0.5;
      \return A reference to the prediction horizon
    */
    double & T() { return T_; }

    /*!
      \return A reference to the descent rate
    */
    double & lam() { return lam_; }

    /*!
      \return A reference to the initial (max) duration for control 
      application used in searching for a valid application duration
    */
    double & maxdt() { return maxdt_; }

    /*!
      \return A reference to the sampling time. This specifies the time 
      between feedback incorporation and the duration between consecutive 
      control calculations
    */
    double & ts() { return ts_; }
    
    /*!
      \return A reference to a 2d vector holding the [max min] saturation 
      values for each element of the control vector
    */
    std::vector<std::vector<double> > & usat() { return usat_; }

    /*!
      \return A reference to the time that each SAC control calculation is 
      assumed to take.  During this interval the previous control is applied 
      while SAC computes the control to apply after this time has elapsed.
    */
    double & calc_tm() { return calc_tm_; }
    
    /*!
      \return A reference to the boolean signifying whether SAC should search
      for an optimal time to apply each control action or if it should apply
      the control associated at the earlies feasible time (at the current time
      plus calc_tm seconds).
    */
    bool & u2search() { return u2search_; }

    /*!
      get using:      mat_type rQ = J_traj.Q();
      get ref using:  mat_type & rQ = J_traj.Q();
      set using:      param.Q() << 1000, 0, 0, 10;
      \return A reference to the mQ_ weight matrix for both getting and 
      setting
    */
    mat_type & Q( ) { return mQ_; }

    /*!
      \return A reference to the mP1_ weight matrix for both getting and 
      setting
    */
    mat_type & P( ) { return mP1_; }

    /*!
      \return A reference to the mR_ weight matrix for both getting and setting
    */
    mat_type & R( ) { return mR_; }

    /*!
      \return A reference to the backtrack_its_ weight matrix for both 
      getting and setting
    */
    size_t & backtrack_its( ) { return backtrack_its_; }

    /*!
      \return A reference to the eps_cost_ for both getting and setting
    */
    double & eps_cost( ) { return eps_cost_; }

    /*!
      \return A reference to mxdes_tf_ for both getting and setting
    */
    vec_type & mxdes_tf( ) { return mxdes_tf_; }
    
  };
  //]

}

#endif // PARAMS_HPP
