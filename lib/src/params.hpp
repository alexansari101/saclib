#ifndef PARAMS_HPP
#define PARAMS_HPP

namespace sac {

  /*!
    User interface for specifying parameters required by SAC 
  */
  class Params {
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
    Eigen::MatrixXd mQ_;  //! incremental quadratic state cost weights 
    Eigen::MatrixXd mP1_; //! terminal quadratic cost weights
    Eigen::MatrixXd mR_;  //! incremental quadratic control cost weights
  public:
    /*!
      Constructor for Params class to hold SAC parameters.
      \param[in] xlen state vector dimension
      \param[in] ulen control vector dimension
    */
    Params( const size_t xlen, 
	    const size_t ulen ) : xlen_(xlen), ulen_(ulen), 
				  T_(.5), lam_(-5), maxdt_(0.2), 
				  ts_(0.01), 
				  usat_(ulen, std::vector<double>(2)), 
				  calc_tm_(ts_), u2search_(false),
				  mQ_(xlen,xlen), mP1_(xlen,xlen),
				  mR_(ulen,ulen)
    { 
      for ( size_t i=0; i<ulen_; i++ ) {
      	usat_[i] = {100, -100}; // default saturation limits
      }
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
      get using:      Eigen::MatrixXd rQ = J_traj.Q();
      get ref using:  Eigen::MatrixXd & rQ = J_traj.Q();
      set using:      param.Q() << 1000, 0, 0, 10;
      \return A reference to the mQ_ weight matrix for both getting and setting
    */
    Eigen::MatrixXd & Q( ) { return mQ_; }

    /*!
      \return A reference to the mP1_ weight matrix for both getting and setting
    */
    Eigen::MatrixXd & P( ) { return mP1_; }

    /*!
      \return A reference to the mR_ weight matrix for both getting and setting
    */
    Eigen::MatrixXd & R( ) { return mR_; }
    
  };
  //]

}

#endif // PARAMS_HPP
