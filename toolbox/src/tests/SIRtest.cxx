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

class SIRtest : public nonlinearSystem {
public:

  SIRtest( int_type neq_in )
  : nonlinearSystem(
      "Semi-implicit approach Example 5",
      "@article{Scheffel:2009,\n"
      "  author  = {Jan Scheffel and Cristian Håkansson},\n"
      "  title   = {Solution of systems of nonlinear equations\n"
      "              – a semi-implicit approach},\n"
      "  journal = {Applied Numerical Mathematics},\n"
      "  volume  = {59},\n"
      "  number  = {10},\n"
      "  pages   = {2430--2443},\n"
      "  year    = {2009},\n"
      "  doi     = {10.1016/j.apnum.2009.05.002},\n"
      "}\n",
      neq_in
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    if ( k == n-1 ) return x(n-1) - 3*cos(x(0));
    return x(k) - cos(x(k+1));
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n-1; ++i )
      f(i) = x(i) - cos(x(i+1));
    f(n-1) = x(n-1) - 3*cos(x(0));
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 2*n;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n-1; ++i ) {
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    SETIJ(n-1,n-1);
    SETIJ(n-1,0);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n-1; ++i ) {
      jac(kk++) = 1;
      jac(kk++) = sin(x(i+1));
    }
    jac(kk++) = 1;
    jac(kk++) = 3*sin(x(0));
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
    if ( n == 2 ) x.fill(-2.0);
    else          x.fill(3.0);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( int_type i = 0; i < n; ++i )
      NONLIN_ASSERT( x(i) > -5 && x(i) < 5, "Bad range" );
  }

  virtual
  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(5);
    L.fill(-5);
  }

};
