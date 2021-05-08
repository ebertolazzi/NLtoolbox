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

class Box3 : public nonlinearSystem {
public:

  Box3()
  : nonlinearSystem(
      "Box3",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      3
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
    f(0) = f(1) = f(2) = 0;
    for ( int_type i = 0; i < 10; ++i ) {
      real_type c = -(i+1) / 10.0;
      real_type fi = exp(c*x(0)) - exp(c*x(1)) - x(2)*(exp(c)-exp(10*c));

      real_type dfidx1 =   c * exp( c * x(0) );
      real_type dfidx2 = - c * exp( c * x(1) );
      real_type dfidx3 = - ( exp(c) - exp( 10*c ) );

      f(0) += 2.0 * fi * dfidx1;
      f(1) += 2.0 * fi * dfidx2;
      f(2) += 2.0 * fi * dfidx3;
    }
  }

  virtual
  int_type
  jacobianNnz() const override
  { return n*n; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();

    for ( int_type i = 0; i < 10; ++i ) {
      real_type c = -(i+1) / 10.0;

      real_type fi = exp(c*x(0)) - exp(c*x(1)) - x(2)*(exp(c) - exp(10*c));

      real_type dfidx1   =   c     * exp( c * x(0) );
      real_type d2fidx11 =   c * c * exp( c * x(0) );
      real_type dfidx2   = - c     * exp( c * x(1) );
      real_type d2fidx22 = - c * c * exp( c * x(1) );
      real_type dfidx3   = - ( exp(c) - exp( 10*c ) );

      jac(faddr(1,1)) += 2.0 * dfidx1 * dfidx1 + 2.0 * fi * d2fidx11;
      jac(faddr(1,2)) += 2.0 * dfidx1 * dfidx2;
      jac(faddr(1,3)) += 2.0 * dfidx1 * dfidx3;

      jac(faddr(2,1)) += 2.0 * dfidx2 * dfidx1;
      jac(faddr(2,2)) += 2.0 * dfidx2 * dfidx2 + 2.0 * fi * d2fidx22;
      jac(faddr(2,3)) += 2.0 * dfidx2 * dfidx3;

      jac(faddr(3,1)) += 2.0 * dfidx3 * dfidx1;
      jac(faddr(3,2)) += 2.0 * dfidx3 * dfidx2;
      jac(faddr(3,3)) += 2.0 * dfidx3 * dfidx3;

    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    x(0) = 1;
    x(1) = 10;
    x(2) = 1;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type  ) const override {
    x(0) = 0;
    x(1) = 10;
    x(2) = 5;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
  }

};
