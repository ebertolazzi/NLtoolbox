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

class DixonFunction : public nonlinearSystem {
public:

  DixonFunction( int_type neq )
  : nonlinearSystem(
      "Dixon Function",
      "@Article{Dixon1988,\n"
      "  author  = {Dixon, L. C. W. and Price, R. C.},\n"
      "  title   = {Numerical experience with the truncated Newton\n"
      "             method for unconstrained optimization},\n"
      "  journal = {Journal of Optimization Theory and Applications},\n"
      "  year    = {1988},\n"
      "  volume  = {56},\n"
      "  number  = {2},\n"
      "  pages   = {245--255},\n"
      "  doi     = {10.1007/BF00939410}\n"
      "}\n",
      neq
    )
  { checkMinEquations(n,2); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    if      ( i == 0   ) return 2*(x(0)-1);
    else if ( i == n-1 ) return 8*n*(2*x(n-1)*x(n-1)-x(n-2))*x(n-1);
    else                 return 8*(i+1)*(2*x(i)*x(i)-x(i-1))-2*(i+2)*(2*x(i+1)*x(i+1)-x(i));
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 2*(x(0)-1);
    for ( int_type i = 1; i < n-1; ++i )
      f(i) = 8*(i+1)*(2*x(i)*x(i)-x(i-1))-2*(i+2)*(2*x(i+1)*x(i+1)-x(i));
    f(n-1) = 8*n*(2*x(n-1)*x(n-1)-x(n-2))*x(n-1);
  }

  virtual
  int_type
  jacobianNnz() const override
  { return 3*n-3; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    ii(kk) = 0; jj(kk) = 0; ++kk;
    for ( int_type i = 1; i < n-1; ++i ) {
      ii(kk) = i; jj(kk) = i-1; ++kk;
      ii(kk) = i; jj(kk) = i;   ++kk;
      ii(kk) = i; jj(kk) = i+1; ++kk;
    }
    ii(kk) = n-1; jj(kk) = n-2; ++kk;
    ii(kk) = n-1; jj(kk) = n-1; ++kk;
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = 2;
    for ( int_type i = 1; i < n-1; ++i ) {
      jac(kk++) = -8*(i+1);
      jac(kk++) = 32*(i+1)*x(i)+2*(i+2);
      jac(kk++) = -8*(i+2)*x(i+1);
    }
    jac(kk++) = -8*n*x(n-1);
    jac(kk++) = 32*n*x(n-1)*x(n-1)+8*n*(2*x(n-1)*x(n-1)-x(n-2));
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 1;
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
