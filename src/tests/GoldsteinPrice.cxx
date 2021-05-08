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

class GoldsteinPrice : public nonlinearSystem {
public:

  GoldsteinPrice()
  : nonlinearSystem(
      "Goldstein Price Polynomial",
      "@book{Michalewicz:1996,\n"
      "  author = {Michalewicz, Zbigniew},\n"
      "  title = {Genetic Algorithms + Data Structures = Evolution Programs (3rd Ed.)},\n"
      "  year = {1996},\n"
      "  isbn = {3-540-60676-9},\n"
      "  publisher = {Springer-Verlag},\n"
      "  address = {Berlin, Heidelberg},\n"
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

    real_type a = x1 + x2 + 1.0;
    real_type b = 19.0 - 14.0 * x1 
                + 3.0 * x1*x1 - 14.0 * x2 
                + 6.0 * x1*x2 + 3.0 * x2*x2;
    real_type c = 2.0 * x1 - 3.0 * x2;
    real_type d = 18.0 - 32.0 * x1 + 12.0 * x1*x1 + 48.0 * x2
                - 36.0 * x1*x2 + 27.0 * x2*x2;

    real_type dbdx1 = -14.0 + 6.0 * (x1+x2);
    real_type dbdx2 = -14.0 + 6.0 * (x1+x2);

    real_type dddx1 = -32.0 + 24.0 * x1 - 36.0 * x2;
    real_type dddx2 = 48.0 - 36.0 * x1 + 54.0 * x2;

    real_type a2 = a*a;
    real_type c2 = c*c;

    f(0) = ( 1.0 + a2 * b ) * (  4.0*c*d + c2 * dddx1 ) 
         + ( 2.0*a*b + a2 * dbdx1 ) * ( 30.0 + c2 * d );
    f(1) = ( 1.0 + a2 * b ) * ( -6.0*c*d + c2 * dddx2 ) 
         + ( 2.0*a*b + a2 * dbdx2 ) * ( 30.0 + c2 * d );
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
    real_type x1 = x(0);
    real_type x2 = x(1);

    real_type a = x1 + x2 + 1.0;

    real_type b = 19.0 - 14.0 * x1 + 3.0 * x1*x1 - 14.0 * x2 + 6.0 * x1 * x2 + 3.0 * x2*x2;
    real_type e = - 14.0 + 6.0 * (x1+x2);
    real_type c = 2.0 * x1 - 3.0 * x2;
    real_type d = 18.0 - 32.0 * x1 + 12.0 * x1*x1 + 48.0 * x2 - 36.0 * x1 * x2 + 27.0 * x2*x2;
    real_type r = - 32.0 + 24.0 * x1 - 36.0 * x2;
    real_type s = 48.0 - 36.0 * x1 + 54.0 * x2;

    real_type a2 = a*a;
    real_type c2 = c*c;

    jac(0) = 2.0 * ( 2.0 * a * b + a2 * e )
           * ( 4.0 * c * d + c2 * r ) + ( 1.0 + a2 * b )
           * ( 8.0 * d + 4.0 * c * r + 4.0 * c * r + 24.0 * c2 )
           + ( 2.0 * b + 2.0 * a * e + 2.0 * a * e + 6.0 * a2 )
           * ( 30.0 + c2 * d );

    jac(1) = ( 2.0 * a * b + a2 * e )
           * ( -2.0 * c * d + c2 * ( r + s ) )
           + ( 1.0 + a2 * b )
           * ( -12.0 * d + 4.0 * c * s -6.0 * c * r - 36.0 * c2 )
           + ( 2.0 * b + 4.0 * a * e + 6.0 * a2 )
           * ( 30.0 + c2 * d );

    jac(2) = ( 2.0 * a * b + a2 * e )
           * ( -2.0 * c * d + c2 * ( s + r ) )
           + ( 1.0 + a2 * b )
           * ( -12.0 * d - 6.0 * c * r + 4.0 * c * s - 36.0 * c2 )
           + ( 2.0 * b + 4.0 * a * e + 6.0 * a2 )
           * ( 30.0 + c2 * d );

    jac(3) = 2.0 * ( 2.0 * a * b + a2 * e )
           * ( -6.0 * c * d + c2 * s )
           + ( 1.0 + a2 * b )
           * ( 18.0 * d - 6.0 * c * s - 6.0 * c * s + 54.0 * c2 )
           + ( 2.0 * b + 2.0 * a * e + 2.0 * a * e + 6.0 * a2 )
           * ( 30.0 + c2 * d );
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = -1;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -0.5;
    x(1) = 0.25;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( std::abs(x(i)) < 10000, "x[" << i << "] = "<< x(i) << " too big");
  }

};
