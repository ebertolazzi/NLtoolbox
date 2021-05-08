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

class SampleProblem18 : public nonlinearSystem {
  real_type const epsi;
public:
  
  SampleProblem18()
  : nonlinearSystem(
     "Sample problem 18",
     "no doc",
     2
    )
  , epsi(1e-8)
  {}
  
  real_type
  fun( real_type x ) const {
    real_type x2 = x*x;
    if ( std::abs(x2) > epsi ) return (1-exp(-x2))/x;
    else                       return x*(1+x2*(x2/6.0-0.5));
  }

  real_type
  fun_1( real_type x ) const {
    real_type x2    = x*x;
    real_type expx2 = exp(-x2);
    if ( std::abs(x2) > epsi ) return 2*expx2-(1-expx2)/x2;
    else                       return 1+x2*((5.0/6.0)*x2-1.5);
  }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return power2(x(1))*fun(x(0));
      case 1: return x(0)*fun(x(1));
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(x(1))*fun(x(0));
    f(1) = x(0)*fun(x(1));
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
    jac(0) = power2(x(1))*fun_1(x(0));
    jac(1) = 2*x(1)*fun(x(0));
    jac(2) = fun(x(1));
    jac(3) = x(0)*fun_1(x(1));
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 2;
    x(1) = 2;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( int_type i = 0; i < n; ++i )
      NONLIN_ASSERT( std::abs(x(i)) < 4, "Bad range" );
  }

  virtual
  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(4);
    L.fill(-4);
  }

};
