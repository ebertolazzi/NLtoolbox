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

class Beale : public nonlinearSystem {
public:

  Beale()
  : nonlinearSystem(
      "Beale",
      "@book{beale1958,\n"
      "  title    = {On an Iterative Method for Finding a Local Minimum\n"
      "              of a Function of More Than One Variable},\n"
      "  author    = {Beale, E.M.L.},\n"
      "  series    = {Technical report\n"
      "               (Princeton University. Statistical Techniques Research Group)},\n"
      "  year      = {1958},\n"
      "  publisher = {Statistical Techniques Research Group,\n"
      "               Section of Mathematical Statistics,\n"
      "               Department of Mathematics, Princeton University}\n"
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
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);

    real_type f1 = 1.5   - x1 * ( 1.0 - x2 );
    real_type f2 = 2.25  - x1 * ( 1.0 - x2 * x2 );
    real_type f3 = 2.625 - x1 * ( 1.0 - x2 * x2 * x2 );

    real_type df1dx1 = x2 - 1;
    real_type df1dx2 = x1;
    real_type df2dx1 = x2 * x2 - 1;
    real_type df2dx2 = 2.0 * x1 * x2;
    real_type df3dx1 = x2 * x2 * x2 - 1.0;
    real_type df3dx2 = 3.0 * x1 * x2 * x2;

    f(0) = 2.0 * ( f1 * df1dx1 + f2 * df2dx1 + f3 * df3dx1 );
    f(1) = 2.0 * ( f1 * df1dx2 + f2 * df2dx2 + f3 * df3dx2 );
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
  jacobian( dvec_t const & x, dvec_t & vals ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);

    real_type f1 = 1.5   - x1 * ( 1.0 - x2 );
    real_type f2 = 2.25  - x1 * ( 1.0 - x2 * x2 );
    real_type f3 = 2.625 - x1 * ( 1.0 - x2 * x2 * x2 );

    real_type df1dx1 = x2 - 1;
    real_type df1dx2 = x1;
    real_type df2dx1 = x2 * x2 - 1;
    real_type df2dx2 = 2.0 * x1 * x2;
    real_type df3dx1 = x2 * x2 * x2 - 1.0;
    real_type df3dx2 = 3.0 * x1 * x2 * x2;

    real_type d2f1dx12 = 1.0;
    real_type d2f1dx21 = 1.0;

    real_type d2f2dx12 = 2.0 * x2;
    real_type d2f2dx21 = 2.0 * x2;
    real_type d2f2dx22 = 2.0 * x1;

    real_type d2f3dx12 = 3.0 * x2 * x2;
    real_type d2f3dx21 = 3.0 * x2 * x2;
    real_type d2f3dx22 = 6.0 * x1 * x2;

    vals(0) = 2.0 * ( df1dx1 * df1dx1 +
                      df2dx1 * df2dx1 +
                      df3dx1 * df3dx1 );

    vals(1) = 2.0 * ( df1dx2 * df1dx1 + f1 * d2f1dx12
                    + df2dx2 * df2dx1 + f2 * d2f2dx12
                    + df3dx2 * df3dx1 + f3 * d2f3dx12 );

    vals(2) = 2.0 * ( df1dx1 * df1dx2 + f1 * d2f1dx21
                    + df2dx1 * df2dx2 + f2 * d2f2dx21
                    + df3dx1 * df3dx2 + f3 * d2f3dx21 );

    vals(3) = 2.0 * ( df1dx2 * df1dx2
                    + df2dx2 * df2dx2 + f2 * d2f2dx22
                    + df3dx2 * df3dx2 + f3 * d2f3dx22 );
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    x(0) = 3;
    x(1) = 0.5;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 1;
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
