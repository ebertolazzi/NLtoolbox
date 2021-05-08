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

class GeneralizedRosenbrock : public nonlinearSystem {
  real_type N;
public:
  
  GeneralizedRosenbrock( int_type neq )
  : nonlinearSystem(
      "Generalized Rosenbrock function",
      "@article{Rosenbrock:1960,\n"
      "  author  = {Rosenbrock, H. H.},\n"
      "  title   = {An Automatic Method for Finding the Greatest\n"
      "             or Least Value of a Function},\n"
      "  journal = {The Computer Journal},\n"
      "  year    = {1960},\n"
      "  volume  = {3},\n"
      "  number  = {3},\n"
      "  pages   = {175--184},\n"
      "  doi = {10.1093/comjnl/3.3.175},\n"
      "}\n\n"
      "@article{More:1981,\n"
      "  author = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  volume = {7},\n"
      "  number = {1},\n"
      "  month = mar,\n"
      "  year = {1981},\n"
      "  pages = {17--41},\n"
      "  doi = {10.1145/355934.355936},\n"
      "}\n",
      neq
    )
  , N(100)
  {
    checkEven(n,2);
  }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type f = 0;
    if ( i == 0 ) {
      f = -4*N*(x(1)-power2(x(0)))*x(0)+2*x(0)-2;
    } else if ( i == n-1 ) {
      f = 2*N*(x(n-1)-power2(x(n-2)));
    } else {
      f = 2*N*(x(i)-power2(x(i-1)))-4*N*(x(i+1)-power2(x(i)))*x(i)+2*x(i)-2;
    }
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = -4*N*(x(1)-power2(x(0)))*x(0)+2*x(0)-2;
    for ( int_type i = 1; i < n-1; ++i )
      f(i) = 2*N*(x(i)-power2(x(i-1)))-4*N*(x(i+1)-power2(x(i)))*x(i)+2*x(i)-2;
    f(n-1) = 2*N*(x(n-1)-power2(x(n-2)));
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    for ( int_type i = 1; i < n-1; ++i ) {
      SETIJ(i,i-1);
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    SETIJ(n-1,n-2);
    SETIJ(n-1,n-1);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = 8*N*power2(x(0))-4*N*(x(1)-power2(x(0)))+2;
    jac(kk++) = -4*N*x(0);
    for ( int_type i = 1; i < n-1; ++i ) {
      jac(kk++) = -4*N*x(i-1);
      jac(kk++) =  2+(12*x(i)*x(i)-4*x(i+1)+2)*N;
      jac(kk++) = -4*N*x(i);
    }
    jac(kk++) = -4*N*x(n-2);
    jac(kk++) = 2*N;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.fill(1);
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; i += 2 ) {
      x(i)   = -1.2;
      x(i+1) = 1;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( std::abs(x(i)) < 10, "Bad range" );
  }

};
