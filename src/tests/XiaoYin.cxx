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

class XiaoYin1 : public nonlinearSystem {
public:

  XiaoYin1()
  : nonlinearSystem(
      "XiaoYin example 3",
      "@article{XiaoYin:2015,\n"
      "  author    = {Xiaoyong Xiao and Hongwei Yin},\n"
      "  title     = {A new class of methods with higher order of convergence for solving systems of nonlinear equations},\n"
      "  Journal   = {Applied Mathematics and Computation},\n"
      "  Number    = {264},\n"
      "  Pages     = {300–-309},\n"
      "  Publisher = {Elsevier},\n"
      "  Year      = {2015},\n"
      "  Doi       = {10.1016/j.amc.2015.04.094},\n"
      "}\n",
      8
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    return x(k)*log(1+x((k+1)%n))-1;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n; ++i )
      f(i) = x(i)*log(1+x((i+1)%n))-1;
  }

  integer
  jacobianNnz() const override
  { return 2*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      ii(kk) = i; jj(kk) = i;       ++kk;
      ii(kk) = i; jj(kk) = (i+1)%n; ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      jac(kk) = log(1+x((i+1)%n));   ++kk;
      jac(kk) = x(i)/(1+x((i+1)%n)); ++kk;
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer  ) const override {
    x.fill(1.5);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class XiaoYin2 : public nonlinearSystem {
public:

  XiaoYin2()
  : nonlinearSystem(
      "XiaoYin example 4",
      "@article{XiaoYin:2015,\n"
      "  author    = {Xiaoyong Xiao and Hongwei Yin},\n"
      "  title     = {A new class of methods with higher order of convergence for solving systems of nonlinear equations},\n"
      "  Journal   = {Applied Mathematics and Computation},\n"
      "  Number    = {264},\n"
      "  Pages     = {300–-309},\n"
      "  Publisher = {Elsevier},\n"
      "  Year      = {2015},\n"
      "  Doi       = {10.1016/j.amc.2015.04.094},\n"
      "}\n",
      16
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    return x(k)*sin(x((k+1)%n))-1;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n; ++i )
      f(i) = x(i)*sin(x((i+1)%n))-1;
  }

  integer
  jacobianNnz() const override
  { return 2*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      ii(kk) = i; jj(kk) = i;       ++kk;
      ii(kk) = i; jj(kk) = (i+1)%n; ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      jac(kk) = sin(x((i+1)%n));      ++kk;
      jac(kk) = x(i)*cos(x((i+1)%n)); ++kk;
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer  ) const override {
    x.fill(-0.85);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class XiaoYin3 : public nonlinearSystem {
public:

  XiaoYin3()
  : nonlinearSystem(
      "XiaoYin example 5",
      "@article{XiaoYin:2015,\n"
      "  author    = {Xiaoyong Xiao and Hongwei Yin},\n"
      "  title     = {A new class of methods with higher order of convergence for solving systems of nonlinear equations},\n"
      "  Journal   = {Applied Mathematics and Computation},\n"
      "  Number    = {264},\n"
      "  Pages     = {300–-309},\n"
      "  Publisher = {Elsevier},\n"
      "  Year      = {2015},\n"
      "  Doi       = {10.1016/j.amc.2015.04.094},\n"
      "}\n",
      35
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    real_type x1  = x(k);
    real_type xi1 = x((k+1)%n);
    return x1*xi1-exp(-x1)-exp(xi1);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n; ++i ) {
      real_type x1  = x(i);
      real_type xi1 = x((i+1)%n);
      f(i) = x1*xi1-exp(-x1)-exp(-xi1);
    }
  }

  integer
  jacobianNnz() const override
  { return 2*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      ii(kk) = i; jj(kk) = i;       ++kk;
      ii(kk) = i; jj(kk) = (i+1)%n; ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      real_type x1  = x(i);
      real_type xi1 = x((i+1)%n);
      jac(kk) = xi1+exp(-x1); ++kk;
      jac(kk) = x1+exp(-xi1); ++kk;
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer  ) const override {
    x.fill(1.2);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
  }

};
