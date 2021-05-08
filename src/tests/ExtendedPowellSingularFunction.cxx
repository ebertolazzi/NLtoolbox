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

class ExtendedPowellSingularFunction : public nonlinearSystem {

  real_type const sqrt5;
  real_type const sqrt10;

public:

  ExtendedPowellSingularFunction()
  : nonlinearSystem(
      "Extended Powell singular function",
      "@article{Powell:1962,\n"
      "  author  = {Powell, M. J. D.},\n"
      "  title   = {An Iterative Method for Finding Stationary\n"
      "             Values of a Function of Several Variables},\n"
      "  journal = {The Computer Journal},\n"
      "  year    = {1962},\n"
      "  volume  = {5},\n"
      "  number  = {2},\n"
      "  pages   = {147--151},\n"
      "  doi     = {10.1093/comjnl/5.2.147}\n"
      "}\n\n"
      "@article{More:1981,\n"
      "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title   = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1981},\n"
      "  volume  = {7},\n"
      "  number  = {1},\n"
      "  pages   = {17--41},\n"
      "  doi     = {10.1145/355934.355936},\n"
      "}\n",
      4
    )
  , sqrt5(sqrt(5.0))
  , sqrt10(sqrt(10.0))
  { checkFour(n,4); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    switch ( i % 4 ) {
      case 0: return x(i+0) + 10 * x(i+1);
      case 1: return sqrt5 * ( x(i+2) - x(i+3) );
      case 2: return power2( x(i+1) - 2 * x(i+2) );
      case 3: return sqrt10 * power2( x(i+0) - x(i+3) );
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n; i += 4 ) {
      f(i+0) = x(i+0) + 10 * x(i+1);
      f(i+1) = sqrt5 * ( x(i+2) - x(i+3) );
      f(i+2) = power2( x(i+1) - 2 * x(i+2) );
      f(i+3) = sqrt10 * power2( x(i+0) - x(i+3) );
    }
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 4*n;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n; i += 4 ) {
      SETIJ(i+0,i+0);
      SETIJ(i+0,i+1);
      SETIJ(i+0,i+2);
      SETIJ(i+0,i+3);

      SETIJ(i+1,i+0);
      SETIJ(i+1,i+1);
      SETIJ(i+1,i+2);
      SETIJ(i+1,i+3);

      SETIJ(i+2,i+0);
      SETIJ(i+2,i+1);
      SETIJ(i+2,i+2);
      SETIJ(i+2,i+3);

      SETIJ(i+3,i+0);
      SETIJ(i+3,i+1);
      SETIJ(i+3,i+2);
      SETIJ(i+3,i+3);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; i += 4 ) {
      jac(kk++) = 1;
      jac(kk++) = 10;
      jac(kk++) = 0;
      jac(kk++) = 0;

      jac(kk++) = 0;
      jac(kk++) = 0;
      jac(kk++) = sqrt5;
      jac(kk++) = -sqrt5;

      jac(kk++) = 0;
      jac(kk++) = 2 * ( x(i+1) - 2 * x(i+2) );
      jac(kk++) = - 4 * ( x(i+1) - 2 * x(i+2) );
      jac(kk++) = 0;

      jac(kk++) = 2 * sqrt10 * ( x(i+0) - x(i+3) );
      jac(kk++) = 0;
      jac(kk++) = 0;
      jac(kk++) = - 2 * sqrt10 * ( x(i+0) - x(i+3) );
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.setZero();
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; i += 4 ) {
      x(i+0) = 3;
      x(i+1) = -1;
      x(i+2) = 0;
      x(i+3) = 1;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
