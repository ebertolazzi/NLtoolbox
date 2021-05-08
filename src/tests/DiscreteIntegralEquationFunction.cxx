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

class DiscreteIntegralEquationFunction : public nonlinearSystem {
public:
  
  DiscreteIntegralEquationFunction( int_type neq )
  : nonlinearSystem(
      "Discrete integral equation function",
      "@article{More:1981,\n"
      "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title   = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1981},\n"
      "  volume  = {7},\n"
      "  number  = {1},\n"
      "  pages   = {17--41},\n"
      "  doi     = {10.1145/355934.355936},\n"
      "}\n",
      neq
    )
  {
    checkMinEquations(n,2);
  }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type h = 1 / real_type ( n + 1 );
    
    real_type tk = real_type ( k + 1 ) / real_type ( n + 1 );
    real_type sum1 = 0;
    for ( int_type j = 0; j < k; ++j ) {
      real_type tj = (j+1) * h;
      sum1 += tj * power3( x(j) + tj + 1 );
    }
    real_type sum2 = 0;
    for ( int_type j = k; j < n; ++j ) {
      real_type tj = (j+1) * h;
      sum2 += (1-tj) * power3( x(j) + tj + 1 );
    }
    return x(k) + h * ( ( 1 - tk ) * sum1 + tk * sum2 ) / 2;
  }
  
  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type h = 1 / real_type ( n + 1 );
    
    for ( int_type k = 0; k < n; ++k ) {
      real_type tk = real_type ( k + 1 ) / real_type ( n + 1 );
      real_type sum1 = 0;
      for ( int_type j = 0; j < k; ++j ) {
        real_type tj = (j+1) * h;
        sum1 += tj * power3( x(j) + tj + 1 );
      }
      real_type sum2 = 0;
      for ( int_type j = k; j < n; ++j ) {
        real_type tj = (j+1) * h;
        sum2 += (1-tj) * power3( x(j) + tj + 1 );
      }
      f(k) = x(k) + h * ( ( 1 - tk ) * sum1 + tk * sum2 ) / 2;
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
    int_type kk = 0;
    for ( int_type k = 0; k < n; ++k ) {
      real_type tk = real_type(k+1) / real_type(n+1);
      for ( int_type j = 0; j < n; ++j ) {
        real_type tj = real_type(j+1) / real_type(n+1);
        real_type temp1 = power2( x(j) + tj + 1 );
        real_type temp2 = min(tk, tj) - tj * tk;
        jac(kk) = 1.5 * temp2 * temp1 / real_type(n+1);
        if ( j == k ) jac(kk) += 1;
        ++kk;
      }
    }
  }

  virtual
  int_type
  numExactSolution() const override {
    if ( n == 2 || n == 5 ) return 1;
    return 0;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    if ( n == 2 ) {
      x(0) = -0.0739748874643476;
      x(1) = -0.162925156390870;
    } else if ( n == 5 ) {
      x(0) = -1.75616539435172;
      x(1) = -5.29228466926422;
      x(2) = -9.72022156411245;
      x(3) = -0.13176631159597;
      x(4) = -0.123717020772933;
    }
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k )
      x(k) = real_type ( (k+1) * (k-n) ) / power2(n+1.0);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
