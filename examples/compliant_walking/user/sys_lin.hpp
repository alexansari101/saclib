#ifndef SYS_LIN_HPP
#define SYS_LIN_HPP

namespace sac {

  //[ Linearizations of the system defined as a class
  // USER SPECIFIED:
  class sys_lin {
    double Cl_;
    double Cr_;
    mat_type DCl_;
    mat_type DCr_;
    int loc_;
    // state_type DzGrndl_;
    // state_type DzGrndr_;
    // double zGrndl_;
    // double zGrndr_;    
    
  public:
    sys_lin( ) : Cl_(0) , Cr_(0) , DCl_(mat_type::Zero(1,8)) , 
		 DCr_(mat_type::Zero(1,8)), loc_(0)
		 /*, DzGrndl_(8) , DzGrndr_(8) , zGrndl_(0) , zGrndr_(0)*/ { }

    void A( const state_type & x, const state_type & u, 
	    mat_type & Amat ) {
      Amat.resize(x.size(),x.size());
      Amat.setZero();
      loc_ = (int) (x[7]+0.5); // pos. floating point location to nearest int
      //
      switch ( loc_ ) {
      case 0:
	Al(x, u, Amat); // Left leg single support dynamics
	break;
      case 1:
	Ad(x, u, Amat); // Double support dynamics
	break;
      case 2:
	Ar(x, u, Amat); // Right leg single support dynamics
	break;
      case 3:
      	Af(x, u, Amat); // flight dynamics
      	break;
      default:
	std::cout << "Location #" << loc_ << " has not been defined\n";
	break;
      }
      //
    }

    void B( const state_type & x, const state_type & u,
	    mat_type & Bmat ) {
      Bmat.resize(x.size(),u.size());
      Bmat.setZero();
      if ( ((int) (x[7]+0.5)) < 3 )
	Bmat(6,0) = 1;
      Bmat(1,1) = 1;
    }
    
    // Left leg single support:
    void Al( const state_type & x, const state_type & /*u*/, 
	     mat_type & Amat ) {
      // zGrndl_ = zGrnd( x[4] );
      // DzGrnd( x[4], DzGrndl_ );
      Cl_ = Cl(x);
      State2Mat( DCl(x), DCl_ );
      //
      Amat(0,1) = 1;
      //
      Amat.row(1) = DCl_ * (x[0]-x[4]) / m ;
      Amat(1,0) += Cl_ / m;
      Amat(1,4) -= Cl_ / m;
      //
      Amat(2,3) = 1;
      //
      Amat.row(3) = DCl_ * (x[2]/*-zGrndl_*/)/m;
      Amat(3,2) += Cl_ / m;
      /* Amat(3,4) -= Cl_ * DzGrndl_ / m */
      //
    }
    // Double support:
    // Right leg single support:
    void Ad( const state_type & x, const state_type & /*u*/, 
	     mat_type & Amat ) {
      // zGrndl_ = zGrnd( x[4] );    zGrndr_ = zGrnd( x[5] );
      // DzGrnd( x[4], DzGrndl_ );   DzGrnd( x[5], DzGrndr_ );
      Cl_ = Cl(x);                Cr_ = Cr(x);
      State2Mat( DCl(x), DCl_ );  State2Mat( DCr(x), DCr_ );
      //
      Amat(0,1) = 1;
      //
      Amat.row(1) = ( DCl_ * (x[0]-x[4]) + DCr_ * (x[0]-x[5]) ) / m ;
      Amat(1,0) += (Cl_ + Cr_) / m;
      Amat(1,4) -= Cl_ / m;
      Amat(1,5) -= Cr_ / m;
      //
      Amat(2,3) = 1;
      //
      Amat.row(3) = ( DCl_ *(x[2]/*-zGrndl_*/) + DCr_ *(x[2]/*-zGrndr_*/) )/ m;
      Amat(3,2) += ( Cl_ + Cr_ ) / m;
      /* Amat(3,4) -= Cl_ * DzGrndl_ / m */
      /* Amat(3,5) -= Cr_ * DzGrndr_ / m */
      //
    }

    // Right leg single support:
    void Ar( const state_type & x, const state_type & /*u*/, 
	     mat_type & Amat ) {
      // zGrndr_ = zGrnd( x[5] );
      // DzGrnd( x[5], DzGrndr_ );
      Cr_ = Cr(x);
      State2Mat( DCr(x), DCr_ );
      //
      Amat(0,1) = 1;
      //
      Amat.row(1) = DCr_ * -(x[5]-x[0]) / m ;
      Amat(1,0) += Cr_ / m;
      Amat(1,5) -= Cr_ / m;
      //
      Amat(2,3) = 1;
      //
      Amat.row(3) = DCr_ * (x[2]/*-zGrndr_*/)/m;
      Amat(3,2) += Cr_ / m;
      /* Amat(3,5) -= Cr_ * DzGrndr_ / m */
      //
    }
    
    // Left leg single support:
    void Af( const state_type & /*x*/, const state_type & /*u*/, 
    	     mat_type & Amat ) {
      Amat(0,1) = 1;
      Amat(2,3) = 1;
    }

  };
  //]


}

#endif  // SYS_LIN_HPP
