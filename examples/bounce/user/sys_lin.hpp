#ifndef SYS_LIN_HPP
#define SYS_LIN_HPP

namespace sac {

  //[ Linearizations of the system defined as a class
  // USER SPECIFIED:
  class sys_lin {

  public:
    sys_lin( ) {  }

    void A( const state_type & /*x*/, const state_type & /*u*/, 
	    mat_type & Amat ) {
      Amat = 0*Amat;
      Amat(0,2) = 1;
      Amat(1,3) = 1; 
    }

    void B( const state_type & /*x*/, const state_type & /*u*/,
	    mat_type & Bmat ) {
      Bmat = Bmat*0;
      Bmat(2,0) = 1;
      Bmat(3,1) = 1;
    }
  };
  //]

}

#endif  // SYS_LIN_HPP
