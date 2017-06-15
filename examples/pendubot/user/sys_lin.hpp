#ifndef SYS_LIN_HPP
#define SYS_LIN_HPP

namespace sac {

  //[ Linearizations of the system defined as a class
  // USER SPECIFIED:
  class sys_lin {
    double g_;
  
  public:
    sys_lin( ) :  g_(9.81) {  }

    void A( const state_type & x, const state_type & u, 
	    mat_type & Amat ) {

      Amat(0,0) = 0; Amat(0,1) = 1; Amat(0,2) = 0; Amat(0,3) = 0;

      Amat(1,0) = (2*pow(p3, 2)*sin(2*(x[0]-x[2])) 
		   * (g_*((-2*p2*p4 + p3*p5)*sin(x[0]) - p3*p5*sin(x[0]-2*x[2])) 
		      - 2*p2*u[0] + 2*p3*sin(x[0]-x[2])*(p3*cos(x[0]-x[2])
							 *pow(x[1], 2) + 
							 p2*pow(x[3], 2))) 
		   + (2*p1*p2 - pow(p3, 2) - pow(p3, 2)*cos(2*(x[0]-x[2])))
		   *(g_*(2*p2*p4 - p3*p5)*cos(x[0]) + g_*p3*p5*cos(x[0]-2*x[2]) 
		     - 2*p3*(p3*cos(2*(x[0]-x[2]))*pow(x[1], 2) 
			     + p2*cos(x[0] - x[2])*pow(x[3], 2)))
		   )/pow( (2*p1*p2 - 2*pow(p3, 2)*pow(cos(x[0] - x[2]), 2)), 2);

      Amat(1,1) = (pow(p3, 2)*sin(2*(x[0]-x[2]))*x[1]
		   )/( -p1*p2 + pow(p3, 2)*pow(cos(x[0]-x[2]), 2) );

      Amat(1,2) = (p3*(g_*(2*pow(p3, 2)*p5*cos(x[0]) 
			   + (2*p2*p3*p4 - 4*p1*p2*p5 + pow(p3, 2)*p5)
			   *cos(x[0]-2*x[2]) 
			   + p3*(-2*p2*p4 + p3*p5)*cos(3*x[0]-2*x[2])) 
		       + 4*p2*p3*sin(2*(x[0]-x[2]))*u[0] - 
		       2*p3*(pow(p3, 2) + (-2*p1*p2 + pow(p3, 2))
			     *cos(2*(x[0]-x[2])))*pow(x[1], 2) 
		       + 2*p2*cos(x[0]-x[2])
		       *(2*p1*p2 - 3*pow(p3, 2) + pow(p3, 2)
			 *cos(2*(x[0]-x[2])))*pow(x[3], 2))
		   )/(4*pow(p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2), 2));

      Amat(1,3) = -((2*p2*p3*sin(x[0]-x[2])*x[3]
		     )/(p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2)));

      Amat(2,0) = 0; Amat(2,1) = 0; Amat(2,2) = 0; Amat(2,3) = 1;

      Amat(3,0) = (p3*(g_*(p3*(p3*p4 - 2*p1*p5)*cos(2.0*x[0]-3.0*x[2]) 
			   + (-4*p1*p2*p4 + pow(p3, 2)*p4 + 2*p1*p3*p5)
			   *cos(2*x[0]-x[2]) + 2*pow(p3, 2)*p4*cos(x[2])) 
		       + 2*(2*p1*p2 + pow(p3, 2) + pow(p3, 2)*cos(2*(x[0]-x[2])))
		       *sin(x[0]-x[2])*u[0] + 2*p1
		       *cos(x[0]-x[2])
		       *(2*p1*p2 - 3*pow(p3, 2) 
			 + pow(p3, 2)*cos(2*(x[0]-x[2])))
		       *pow(x[1], 2) - 2*p3*(pow(p3, 2) + (-2*p1*p2 + pow(p3, 2))
					     *cos(2*(x[0]-x[2])))*pow(x[3], 2))
		   )/(4*pow((p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2) ), 2));

      Amat(3,1) = (2*p1*p3*sin(x[0]-x[2])*x[1]
		   )/(p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2));

      Amat(3,2) = (2*pow(p3, 2)*cos(x[0]-x[2])*sin(x[0]-x[2])
		   *(-g_*p3*p4*cos(x[0]-x[2])*sin(x[0]) 
		     + g_*p1*p5*sin(x[2]) - p3*cos(x[0]-x[2])*u[0] 
		     + p3*sin(x[0]-x[2])
		     *(p1*pow(x[1], 2) + p3*cos(x[0]-x[2])*pow(x[3], 2))) 
		   -1/4.0*(2*p1*p2 - pow(p3, 2) - pow(p3, 2)*cos(2*(x[0]-x[2])))
		   *(-g_*p3*p4*cos(2*x[0]-x[2]) + g_*(p3*p4 - 2*p1*p5)
		     *cos(x[2]) + 2*p3*(sin(x[0]-x[2])*u[0] 
					+ p1*cos(x[0]-x[2])*pow(x[1], 2) 
					+ p3*cos(2*(x[0]-x[2]))
					*pow(x[3], 2)))
		   )/pow( (p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2)) , 2);

      Amat(3,3) = -((pow(p3, 2)*sin(2*(x[0]-x[2]))*x[3]
		     )/(-p1*p2 + pow(p3, 2)*pow(cos(x[0]-x[2]), 2)));
    }

    void B( const state_type & x, const state_type & /*u*/, 
	    mat_type & Bmat ) {
      Bmat << 0,   
	p2/( p1*p2 - pow(p3, 2)*pow(cos(x[0]-x[2]), 2) ),
	0,
	1.0/( p3*cos(x[0]-x[2]) - ( p1*p2/cos(x[0]-x[2]) )/p3 );
    }
  };
  //]

}

#endif  // SYS_LIN_HPP
