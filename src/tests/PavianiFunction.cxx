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

class PavianiFunction : public nonlinearSystem {

public:

  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  PavianiFunction()
  : nonlinearSystem(
      "Paviani function",
      "@book{himmelblau:1972,\n"
      "  author    = {Himmelblau, D.M.},\n"
      "  title     = {Applied nonlinear programming},\n"
      "  year      = {1972},\n"
      "  publisher = {McGraw-Hill}\n"
      "}\n",
      10
    )
  {}

  real_type
  loglog( real_type x ) const {
    return power2(log(10-x))+power2(log(x-2));
  }

  real_type
  loglog_D( real_type x ) const {
    real_type t1 = 10 - x;
    real_type t2 = x - 2;
    return 2*(log(t2)/t2-log(t1)/t1);
  }

  real_type
  loglog_DD( real_type x ) const {
    real_type t1 = 10 - x;
    real_type t2 = x - 2;
    return 2*((1-log(t1))/(t1*t1)+(1-log(t2))/(t2*t2));
  }

  real_type
  mul( dvec_t const & x ) const {
    real_type res = 1;
    for ( integer i = 0; i < 10; ++i )
      res *= power2(abs(x(i)));
    return res;
  }

  real_type
  mul_D( dvec_t const & x, integer k ) const {
    real_type res = 1;
    for ( integer i = 0; i < 10; ++i ) {
      res *= x(i);
      if ( i == k ) res *= 2;
      else          res *= x(i);
    }
    return res;
  }

  real_type
  mul_DD( dvec_t const & x, integer i, integer j ) const {
    real_type res = 1;
    if ( i == j ) {
      for ( integer k = 0; k < 10; ++k ) {
        if ( i == k ) res *= 2;
        else          res *= power2(x(k));
      }
    } else {
      for ( integer k = 0; k < 10; ++k ) {
        res *= x(k);
        if ( i == k || j == k ) res *= 2;
        else                    res *= x(k);
      }
    }
    return res;
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    return loglog_D(x(k)) - mul_D(x,k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < 10; ++i ) {
      if ( x(i) > 2 && x(i) < 10 )
        f(i) = loglog_D(x(i)) - mul_D(x,i);
      else
        f(i) = nan("PavianiFunction");
    }
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
    integer kk = 0;
    for ( integer i = 0; i < 10; ++i ) {
      for ( integer j = 0; j < 10; ++j ) {
        jac(kk) = -mul_DD(x,i,j);
        if ( i == j ) jac(kk) += loglog_DD(x(i));
        ++kk;
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    /*x(0) = 2.001;
    x(1) = 9.999;
    x(2) = 2.001;
    x(3) = 9.999;
    x(4) = 2.001;
    x(5) = 9.999;
    x(6) = 2.001;
    x(7) = 9.999;
    x(8) = 2.001;
    x(9) = 9.999;
    */
    x(0) = 2.1;
    x(1) = 9.9;
    x(2) = 2.1;
    x(3) = 9.9;
    x(4) = 2.1;
    x(5) = 9.9;
    x(6) = 2.1;
    x(7) = 9.9;
    x(8) = 2.1;
    x(9) = 9.9;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < 10; ++i ) {
      UTILS_ASSERT(
        x(i) > 2 && x(i) < 10,
        "x[{}] = {} must be in (2,10)", i, x(i)
      );
    }
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(10);
    L.fill(2);
  }

};
