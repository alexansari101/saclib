#ifndef SYS_LIN_HPP
#define SYS_LIN_HPP

namespace sac {

  //[ Linearizations of the system defined as a class
  // USER SPECIFIED:
  class sys_lin {
    state_type Dgrnd_;
    double zGrnd_;
    double L_;
    //
    double dzGrnddxt_;
    double dzGrnddyt_;
    double dLdzGrnd_;
    double dLdx_;
    double dLdy_;
    double dLdz_;
    double dLdxt_;
    double dLdyt_;
    double denom_;
    //
    int loc_;
    //
    
  public:
    sys_lin( ) : Dgrnd_(8) , zGrnd_(0) , L_(0) ,
                 dzGrnddxt_(0) , dzGrnddyt_(0) , dLdzGrnd_(0) ,
                 dLdx_(0), dLdy_(0), dLdz_(0), dLdxt_(0), 
		 dLdyt_(0), denom_(1), loc_(0) {  }

    void A( const state_type & x, const state_type & u, 
	    mat_type & Amat ) {
      Amat.resize(x.size(),x.size());
      Amat.setZero();
      loc_ = (int) (x[8]+0.5); // pos. floating point location to nearest int
      zGrnd_ = zGrndToe( x[6], x[7] );
      L_ = pow( pow(x[0]-x[6],2)+pow(x[2]-x[7],2)
		+pow(x[4]-zGrnd_,2) ,0.5);
      DzGrndToe( x, Dgrnd_ );
      //
      dzGrnddxt_ = Dgrnd_[6];
      dzGrnddyt_ = Dgrnd_[7];
      dLdzGrnd_ = (zGrnd_-x[4])/L_;
      dLdx_ = (x[0]-x[6])/L_;
      dLdy_ = (x[2]-x[7])/L_;
      dLdz_ = (x[4]-zGrnd_)/L_;
      dLdxt_ = -dLdx_ + dLdzGrnd_*dzGrnddxt_;
      dLdyt_ = -dLdy_ + dLdzGrnd_*dzGrnddyt_;
      //
      Amat = 0*Amat;
      if ( loc_ > 0 ) { // flight linearization
	Amat(0,1) = 1;
	Amat(2,3) = 1;
	Amat(4,5) = 1;
	Amat(6,1) = 1;
	Amat(7,3) = 1;
      } else {             // stance linearization
	denom_ = m*pow(L_,2);
	//
	Amat(0,1) = 1;
	//
	Amat(1,0) = (-k*pow(L_,2) + L_*(k*L0 + u[2]) 
		     - (k*L0 + u[2])*(x[0] - x[6])*dLdx_ )
		     / (denom_);
	Amat(1,2) = -(k*L0 + u[2])*(x[0] - x[6])*dLdy_ / (denom_);
	Amat(1,4) = -(k*L0 + u[2])*(x[0] - x[6])*dLdz_ / (denom_);
	Amat(1,6) = (k*pow(L_,2) - L_*(k*L0 + u[2]) 
		     - (k*L0 + u[2])*(x[0] - x[6])*dLdxt_) 
	            / (denom_);
	Amat(1,7) = -(k*L0 + u[2])*(x[0] - x[6])*dLdyt_ / (denom_);
	//
	Amat(2,3) = 1;
	//
	Amat(3,0) = -(k*L0 + u[2])*(x[2] - x[7])*dLdx_ / (denom_);
	Amat(3,2) = (-k*pow(L_,2) + L_*(k*L0 + u[2]) 
		     - (k*L0 + u[2])*(x[2] - x[7])*dLdy_ )
		     / (denom_);
	Amat(3,4) = -(k*L0 + u[2])*(x[2] - x[7])*dLdz_ / (denom_);
	Amat(3,6) = -(k*L0 + u[2])*(x[2] - x[7])*dLdxt_ / (denom_);
	Amat(3,7) = (k*pow(L_,2) - L_*(k*L0 + u[2]) 
		     - (k*L0 + u[2])*(x[2] - x[7])*dLdyt_) 
	            / (denom_);
	//
	Amat(4,5) = 1;
	//
	Amat(5,0) = (k*L0 + u[2])*(zGrnd_ - x[4])*dLdx_ / (denom_);
	Amat(5,2) = (k*L0 + u[2])*(zGrnd_ - x[4])*dLdy_ / (denom_);
	Amat(5,4) = (-k*pow(L_,2) + L_*(k*L0 + u[2]) 
		     + (k*L0 + u[2])*(zGrnd_ - x[4])*dLdz_ )
		     / (denom_);
	Amat(5,6) = (k*pow(L_,2)*dzGrnddxt_ - L_*(k*L0 + u[2])*dzGrnddxt_
		     + (k*L0 + u[2])*(zGrnd_ - x[4])*dLdxt_) 
	            / (denom_);
	Amat(5,7) = (k*pow(L_,2)*dzGrnddyt_ - L_*(k*L0 + u[2])*dzGrnddyt_
		     + (k*L0 + u[2])*(zGrnd_ - x[4])*dLdyt_) 
	            / (denom_);
      }
    }

    void B( const state_type & x, const state_type & u,
	    mat_type & Bmat ) {
      Bmat.resize(x.size(),u.size());
      Bmat.setZero();
      //
      zGrnd_ = zGrndToe( x[6], x[7] );
      L_ = pow( pow(x[0]-x[6],2)+pow(x[2]-x[7],2)
		+pow(x[4]-zGrnd_,2) ,0.5);
      if ( ((int) (x[8]+0.5)) > 0 ) { // flight linearization
	Bmat(6,0) = 1;
	Bmat(7,1) = 1;
      } else {             // stance linearization
	Bmat(1,2) = (x[0]-x[6])/(m*L_);
	Bmat(3,2) = (x[2]-x[7])/(m*L_);
	Bmat(5,2) = (x[4]-zGrnd_)/(m*L_);
      }
    }

  };
  //]


}

#endif  // SYS_LIN_HPP
