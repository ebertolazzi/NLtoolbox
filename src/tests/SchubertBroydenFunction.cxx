/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

/*
Numerical Performance of Abs Codes for Nonlinear Systems of Equations
E. Bodon (1), A. Del Popolo (1, 2, 3), L. Luksan (4), E. Spedicato (1)
2001
https://arxiv.org/abs/math/0106029
*/

class SchubertBroydenFunction : public nonlinearSystem {

public:

  SchubertBroydenFunction( int_type neq )
  : nonlinearSystem(
      "Schubert Broyden function",
      "no doc",
      neq
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    if ( i == 0   ) return (3-x(0))*x(0)+1-2*x(1);
    if ( i == n-1 ) return (3-x(n-1))*x(n-1)+1-2*x(n-2);
    return (3-x(i))*x(i)-x(i-1)-2*x(i+1);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0)   = (3-x(0))*x(0)+1-2*x(1);
    f(n-1) = (3-x(n-1))*x(n-1)+1-2*x(n-2);
    for ( int_type i = 1; i < n-1; ++i )
      f(i) = (3-x(i))*x(i)-x(i-1)-2*x(i+1);
  }

  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(n-1,n-1);
    SETIJ(n-1,n-2);
    for ( int_type i = 1; i < n-1; ++i ) {
      SETIJ(i,i-1);
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;

    jac(kk++) = 3-2*x(0);
    jac(kk++) = -2;
    jac(kk++) = 3-2*x(n-1);
    jac(kk++) = -2;

    for ( int_type i = 1; i < n-1; ++i ) {
      jac(kk++) = -1;
      jac(kk++) = 3-2*x(i);
      jac(kk++) = -2;
    }

  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = -1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
