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

class DeistSefor : public nonlinearSystem {
  real_type beta[6];

public:

  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  DeistSefor()
  : nonlinearSystem(
      "DeistSefor function",
      "@Article{Martinez1980,\n"
      "  author  = {Mart{\\'i}nez, Jos{\\'e} Mario},\n"
      "  title   = {Solving nonlinear simultaneous equations with\n"
      "             a generalization of Brent's method},\n"
      "  journal = {BIT Numerical Mathematics},\n"
      "  year    = {1980},\n"
      "  volume  = {20},\n"
      "  number  = {4},\n"
      "  pages   = {501--510},\n"
      "  doi     = {10.1007/BF01933643}\n"
      "}\n",
      6
    )
  {
    beta[0] = 0.02249;
    beta[1] = 0.02166;
    beta[2] = 0.02083;
    beta[3] = 0.02;
    beta[4] = 0.01918;
    beta[5] = 0.01835;
  }
  
  real_type cot( real_type x ) const { return 1/tan(x); }
  real_type csc2( real_type x ) const { return 1/power2(sin(x)); }

  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type f = 0;
    for ( int_type j = 0; j < n; ++j )
      if ( i != j ) f += cot(beta[i]*x(j));
    return f;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n; ++i ) {
      f(i) = 0;
      for ( int_type j = 0; j < n; ++j )
        if ( i != j ) f(i) += cot(beta[i]*x(j));
    }
  }

  int_type
  jacobianNnz() const override
  { return n*(n-1); }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        if ( i != j )
          { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = 0; j < n; ++j )
        if ( i != j )
          jac(kk++) = -beta[i]/power2(sin(beta[i]*x(j)));
    }
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numInitialPoint() const override { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(75);
  }

};
