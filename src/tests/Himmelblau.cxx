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

class Himmelblau : public nonlinearSystem {
public:

  Himmelblau()
  : nonlinearSystem(
      "Himmelblau function",
      "@book{himmelblau:1972,\n"
      "  author    = {Himmelblau, D.M.},\n"
      "  title     = {Applied nonlinear programming},\n"
      "  year      = {1972},\n"
      "  publisher = {McGraw-Hill}\n"
      "}\n\n"
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

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type x0_2 = x(0)*x(0);
    real_type x1_2 = x(1)*x(1);
    if ( k == 0 )
      return 4 * ( x0_2 + x(1) - 11 ) * x(0) + 2 * ( x(0) + x1_2 - 7 );
    else
      return 2 * ( x0_2 + x(1) - 11 )        + 4 * ( x(0) + x1_2 - 7 ) * x(1);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x0_2 = x(0)*x(0);
    real_type x1_2 = x(1)*x(1);
    f(0) = 4 * ( x0_2 + x(1) - 11 ) * x(0) + 2 * ( x(0) + x1_2 - 7 );
    f(1) = 2 * ( x0_2 + x(1) - 11 )        + 4 * ( x(0) + x1_2 - 7 ) * x(1);
  }

  virtual
  int_type
  jacobianNnz() const override
  { return 4; }

  virtual
  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = 8 * x(0) * x(0) + 4 * ( x(0)*x(0) + x(1) - 11 ) + 2;
    jac(1) = 4 * x(0) + 4 * x(1);
    jac(2) = 4 * x(0) + 4 * x(1);
    jac(3) = 2        + 8 * x(1) * x(1) + 4 * ( x(0) + x(1) * x(1) - 7 );
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 3.0;
    x(1) = 2.0;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -1.3;
    x(1) =  2.7;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
