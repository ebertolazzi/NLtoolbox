
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

class BertolazziRootPlusSquare : public nonlinearSystem {
  real_type const epsilon;
  real_type const delta;
  real_type       xmin;
public:

  BertolazziRootPlusSquare()
  : nonlinearSystem(
      "Bertolazzi: root+square.",
      "no doc",
      2
    )
  , epsilon(1E-6)
  , delta(1E-8)
  {
    xmin = epsilon*(1.0/exp(1.0)-1.0);
  }

  real_type
  evalFk( dvec_t const & x_in, integer k ) const override {
    real_type x = x_in(0);
    real_type y = x_in(1);
    switch ( k ) {
      case 0: return x > xmin ? log1p(log1p(x/epsilon)) : nan("nan");
      case 1: return delta*y+power2(y);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x = x_in(0);
    real_type y = x_in(1);
    f(0) = x > xmin ? log1p(log1p(x/epsilon)) :nan("nan");
    f(1) = delta*y+power2(y);
  }

  integer
  jacobianNnz() const override
  { return 2; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 1; j(1) = 1;
  }

  void
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x = x_in(0);
    real_type y = x_in(1);
    jac(0) = x > xmin ? 1/((x+epsilon)*(1+log1p(x/epsilon))) : nan("nan");
    jac(1) = delta+2*y;
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1000;
    x(1) = 1000;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible ( dvec_t const & x ) const override {
    UTILS_ASSERT(
      x(0) >= xmin,
      "BertolazziRootPlusSquare: x = {} must be >= {}", x(0), xmin
    );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    L(0) = xmin;
    L(1) = -real_max;
    U(0) = U(1) = real_max;
  }

};

/*
 * Reference:
 *  E.Bertolazzi
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BertolazziAtanPlusQuadratic : public nonlinearSystem {
  real_type const epsilon;
public:

  BertolazziAtanPlusQuadratic()
  : nonlinearSystem(
      "Bertolazzi: atan+quadratic",
      "no doc",
      2
    )
  , epsilon(1E-9)
  { }

  real_type
  evalFk( dvec_t const & x_in, integer k ) const override {
    real_type x = x_in(0);
    real_type y = x_in(1);
    switch ( k ) {
     case 0: return atan(x/epsilon);
     case 1: return epsilon*y+x*y;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x = x_in(0);
    real_type y = x_in(1);
    f(0) = atan(x/epsilon);
    f(1) = epsilon*y+x*y;
  }

  integer
  jacobianNnz() const override
  { return 3; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 1; j(1) = 0;
    i(2) = 1; j(2) = 1;
  }

  void
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x  = x_in(0);
    real_type y  = x_in(1);
    real_type xe = x/epsilon;
    jac(0) = 1/(epsilon+x*xe);
    jac(1) = y;
    jac(2) = epsilon+x;
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = 0;
    x(1) = 0;
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
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT(
        abs(x(i)) < 10,
        "x[{}] = {} out of range [-10,10]", i, x(i)
      );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    L.fill(-10);
    U.fill(10);
  }

};

/*
 * Reference:
 *  E.Bertolazzi
 *
 * f = x*(1+x^2)*exp(-x)+exp(x/10)-1
 * g = x*exp(-x)+exp(x/10)-(3*exp(-3)+exp(3/10));
 *
 * f-g
 * f+g
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BertolazziHard : public nonlinearSystem {
public:

  BertolazziHard()
  : nonlinearSystem(
      "Bertolazzi Hard.",
      "no doc",
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x_in, integer k ) const override {
    real_type x = x_in(0);
    real_type y = x_in(1);
    real_type t1 = exp(x/10);
    real_type t2 = exp(y/10);
    real_type t3 = exp(3.0/10.0);
    real_type t4 = x*(x*x+1)*exp(-x);
    real_type t5 = y*exp(-y);
    real_type t6 = 3*exp(-3);
    switch ( k ) {
     case 0: return t6 - 1 + t4 - t5 + t1 - t2 + t3;
     case 1: return -t6 - 1 + t4 + t5 + t1 + t2 - t3;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x = x_in(0);
    real_type y = x_in(1);
    real_type t1 = exp(x/10);
    real_type t2 = exp(y/10);
    real_type t3 = exp(3.0/10.0);
    real_type t4 = x*(x*x+1)*exp(-x);
    real_type t5 = y*exp(-y);
    real_type t6 = 3*exp(-3);
    f(0) = t6 - 1 + t4 - t5 + t1 - t2 + t3;
    f(1) = -t6 - 1 + t4 + t5 + t1 + t2 - t3;
  }

  integer
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
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x  = x_in(0);
    real_type y  = x_in(1);
    real_type t1 = x*x;
              t1 = exp(x/10)/10+exp(-x)*(3*t1+1-x*(1+t1));
    real_type t2 = -exp(y/10)/10 + exp(-y) *(y-1);
    jac(0) = t1;
    jac(1) = t2;
    jac(2) = t1;
    jac(3) = -t2;
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 10;
    x(1) = -1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*
 * Reference:
 *  E.Bertolazzi
 *
 * f = x*(1+x^2)*exp(-x)+exp(x/10)-1
 */

class BertolazziSingleEQ : public nonlinearSystem {
public:

  BertolazziSingleEQ()
  : nonlinearSystem(
      "Bertolazzi Single EQ.",
      "no doc",
      1
    )
  { }

  real_type
  evalFk( dvec_t const & x_in, integer k ) const override {
    real_type x = x_in(0);
    return exp(-x) * (x*x + 1) * x + exp(x/10) - 1;
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x = x_in(0);
    f(0) = exp(-x) * (x*x + 1) * x + exp(x/10) - 1;
  }

  integer
  jacobianNnz() const override
  { return 1; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
  }

  void
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x = x_in(0);
    real_type t1 = x * x;
    real_type t2 = exp(-x);
    real_type t9 = exp(x/10);
    jac(0) = 3*t2*t1+t2-t1*x*t2-t2*x+t9/10;
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = 0;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 10;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
