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

class BrownBadlyScaled : public nonlinearSystem {
public:

  BrownBadlyScaled()
  : nonlinearSystem(
      "Brown Badly Scaled Function",
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
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    
    f(0) = (2 * x1 * ( 1 + x2 * x2 ) - 4 * x2) - 2000000.0;
    f(1) = (2 * x2 * ( 1 + x1 * x1 ) - 4 * x1) - 0.000004;

  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);

    jac(0) = 2*( 1 + x2 * x2 );
    jac(1) = jac(2) = 4 * ( x1 * x2 - 1 );
    jac(3) = 2*( 1 + x1 * x1 );
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 1e6;
    x(1) = 2e-6;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
