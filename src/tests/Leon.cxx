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

class Leon : public nonlinearSystem {
public:

  Leon()
  : nonlinearSystem(
      "Leon cubic valley function",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type x0_3 = x(0)*x(0)*x(0);
    if ( k == 0 )
      return - 600.0 * ( x(1) - x0_3 ) * x(0) * x(0) - 2.0 * ( 1.0 - x(0) );
    else
      return 200.0 * ( x(1) - x0_3 );
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x0_3 = x(0)*x(0)*x(0);
    f(0) = - 600.0 * ( x(1) - x0_3 ) * x(0) * x(0) - 2.0 * ( 1.0 - x(0) );
    f(1) = 200.0 * ( x(1) - x0_3 );
  }

  int_type
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x0_4 = x(0)*x(0)*x(0)*x(0);
    jac(0) = - 1200.0 * x(0) * x(1) + 3000.0 * x0_4 + 2.0;
    jac(1) = - 600.0 * x(0) * x(0);
    jac(2) = - 600.0 * x(0) * x(0);
    jac(3) =   200.0;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1.0;
    x(1) = 1.0;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -1.2;
    x(1) = -1.0;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i ) {
    //  NONLIN_ASSERT( abs(x(i)) < 1000, "X!!!!" );
    //}
  }

};
