#ifndef SETTING_HPP
#define SETTING_HPP

namespace sac {

  /*********************************************/
  /* Parameters */
  double T = 0.28;  // 1.2; 2.0 // prediction horizon
  double lam = -10; 
  double maxdt = .2;
  double ts = 0.001; // .0167;
  double usat[1][2] = { {25, -25} };
  double calc_tm = ts;
  /**/

  const size_t xlen = 4, ulen = 1;
  bool u2Search = false;

}


namespace sac {
  /*!
    User interface for specifying parameters required by SAC 
   */
  class Params {
    const size_t xlen_; //! state dimension
    const size_t ulen_; //! control dimension
    double T_;          //! prediction horizon
    double lam_;        //! descent rate (should be < 0)
    double maxdt_;      //! initial (max) duration for control application
    double ts_;         //! sampling time (duration between control calcs)
    //    double *usat_;      //! control saturation ranges
    double calc_tm_;    //! time assumed for SAC calcs
    bool u2search_;     //! search for best time to apply control?
    /* cost parameters */
    // Eigen::Matrix< double, xlen_, xlen_ > mQ_;
    // Eigen::Matrix< double, xlen_, xlen_ > mP1_;
    // Eigen::Matrix< double, ulen_, ulen_ > mR_;
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
				  // usat_(new double(ulen * 2)), 
				  calc_tm_(ts_), u2search_(false) 
    { 
      // for ( size_t i=0; i<ulen; i++ ) {
      // 	usat_[i][0] = 100;  // default saturation limits
      // 	usat_[i][1] = -100; // default saturation limits
      // }
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
    
    // /*!
    //   \return A reference to the initial (max) duration for control 
    //   application used in searching for a valid application duration
    // */
    // double & usat() { return usat_; }

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

    // /*!
    //   get using:      Eigen::Matrix< double, xlen, xlen > rQ = J_traj.Q();
    //   get ref using:  Eigen::Matrix< double, xlen, xlen > & rQ = J_traj.Q();
    //   set using:      param.Q() << 1000, 0, 0, 10;
    //   \return A reference to the mQ_ weight matrix for both getting and setting
    // */
    // Eigen::Matrix< double, xlen, xlen > & Q( ) { return mQ_; }

    // /*!
    //   \return A reference to the mP1_ weight matrix for both getting and setting
    // */
    // Eigen::Matrix< double, xlen, xlen > & P( ) { return mP1_; }

    // /*!
    //   \return A reference to the mR_ weight matrix for both getting and setting
    // */
    // Eigen::Matrix< double, ulen, ulen > & R( ) { return mR_; }
    
  };
}

#endif  // SETTING_HPP
