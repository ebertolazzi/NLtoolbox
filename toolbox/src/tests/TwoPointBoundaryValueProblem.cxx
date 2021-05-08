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

class TwoPointBoundaryValueProblem : public nonlinearSystem {

public:

  TwoPointBoundaryValueProblem( int_type neq )
  : nonlinearSystem(
      "Two-Point Boundary Value Problem",
      "@article{More:1979,\n"
      "  author  = {Mor{\'e}, Jorge J. and Cosnard, Michel Y.},\n"
      "  title   = {Numerical Solution of Nonlinear Equations},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1979},\n"
      "  volume  = {5},\n"
      "  number  = {1},\n"
      "  pages   = {64--85},\n"
      "  doi     = {10.1145/355815.355820},\n"
      "}\n\n",
      neq
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    if ( i == 0   ) return x(0);
    if ( i == n-1 ) return x(n-1);
    real_type h = 1.0/(n-1.0);
    real_type t = h*i;
    return 2*x(i)-x(i-1)-x(i+1)+(h*h/2)*power3(x(i)+t+1);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type h = 1.0/(n-1.0);
    
    f(0)   = x(0);
    f(n-1) = x(n-1);
    for ( int_type i = 1; i < n-1; ++i ) {
      real_type t = h*i;
      f(i) = 2*x(i)-x(i-1)-x(i+1)+(h*h/2)*power3(x(i)+t+1);
    }
  }

  int_type
  jacobianNnz() const override {
    return 3*(n-2)+2;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(n-1,n-1);
    for ( int_type i = 1; i < n-1; ++i ) {
      SETIJ(i,i-1);
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type  kk = 0;
    real_type h  = 1.0/(n-1.0);
    jac(kk++) = 1;
    jac(kk++) = 1;
    for ( int_type i = 1; i < n-1; ++i ) {
      real_type t = h*i;
      jac(kk++) = -1;
      jac(kk++) = 2 + 1.5*(h*h)*power2(x(i)+t+1);
      jac(kk++) = -1;
    }
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    real_type h = 1.0/(n-1.0);
    for ( int_type i = 0; i < n; ++i ) {
      real_type t = h*i;
      x(i) = t*(t-1);
    }
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
