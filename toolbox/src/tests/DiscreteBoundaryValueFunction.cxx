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

class DiscreteBoundaryValueFunction : public nonlinearSystem {
   real_type h;
public:

  DiscreteBoundaryValueFunction( integer neq )
  : nonlinearSystem(
      "Discrete boundary value function",
      "@article{More:1979,\n"
      "  author  = {Mor{\'e}, Jorge J. and Cosnard, Michel Y.},\n"
      "  title   = {Numerical Solution of Nonlinear Equations},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1979},\n"
      "  volume  = {5},\n"
      "  number  = {1},\n"
      "  pages   = {64--85},\n"
      "  doi     = {10.1145/355815.355820},\n"
      "}\n\n"
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
  , h(1/ real_type(n+1))
  {
    checkMinEquations(n,1);
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    real_type f = 2*x(k) + 0.5 * power2(h) * power3( (x(k)+1) + (k+1) * h );
    if ( k > 0   ) f -= x(k-1);
    if ( k < n-1 ) f -= x(k+1);
    return f;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer k = 0; k < n; ++k ) {
      f(k) = 2*x(k) + 0.5 * power2(h) * power3( (x(k)+1) + (k+1) * h );
      if ( k > 0   ) f(k) -= x(k-1);
      if ( k < n-1 ) f(k) -= x(k+1);
    }
  }

  integer
  jacobianNnz() const override
  { return 3*n-2; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;

    for ( integer i = 0; i < n; ++i )
      { ii(kk) = jj(kk) = i; ++kk; }

    for ( integer i = 0; i < n-1; ++i ) {
      ii(kk) = i+1; jj(kk) = i;   ++kk;
      ii(kk) = i;   jj(kk) = i+1; ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      jac(kk++) = 2 + 1.5*h*h*power2( x(i) + h*(i+1) + 1 );

    for ( integer i = 0; i < n-1; ++i ) {
      jac(kk++) = -1;
      jac(kk++) = -1;
    }
  }

  integer
  numExactSolution() const override {
    if ( n == 2 || n == 5 ) return 1;
    return 0;
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    if ( n == 2 ) {
      x(0) = -0.128246763033732;
      x(1) = -0.159267567244641;
    } else if ( n == 5 ) {
      x(0) = -0.0750221292923205;
      x(1) = -0.131976210352191;
      x(2) = -0.164848771909337;
      x(3) = -0.164664680215801;
      x(4) = -0.117417651684194;
    }
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    for ( integer k = 0; k < n; ++k )
      x(k) = real_type ( (k+1) * (k-n) ) / power2(n+1.0);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
