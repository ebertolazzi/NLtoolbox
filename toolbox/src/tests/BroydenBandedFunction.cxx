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

class BroydenBandedFunction : public nonlinearSystem {
  integer const ml;
  integer const mu;
public:
  
  BroydenBandedFunction()
  : nonlinearSystem(
      "Broyden Banded Function",
      "@article{Broyden:1971,\n"
      "  author  = {Broyden, C. G.},\n"
      "  title   = {The convergence of an algorithm for solving sparse nonlinear systems},\n"
      "  journal = {Mathematics of Computation},\n"
      "  volume  = {25},\n"
      "  year    = {1971},\n"
      "  pages   = {285--294},\n"
      "  doi     = {10.2307/2004922},\n"
      "}\n\n"
      "@article{More:1981,\n"
      "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title   = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  volume  = {7},\n"
      "  number  = {1},\n"
      "  year    = {1981},\n"
      "  pages   = {17--41},\n"
      "  doi     = {10.1145/355934.355936},\n"
      "}\n",
      10
    )
  , ml(5)
  , mu(1)
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
  	integer k1 = max(0,  k-ml);
    integer k2 = min(n-1,k+mu);
    real_type temp = 0;
    for ( integer j = k1; j <= k2; ++j ) {
      if ( j != k ) temp += x(j)*(1+x(j));
    }
    return x(k)*(2+5*x(k)*x(k)) + 1 - temp;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer k = 0; k < n; ++k ) {
    	integer k1 = max(0,  k-ml);
      integer k2 = min(n-1,k+mu);
      real_type temp = 0;
      for ( integer j = k1; j <= k2; ++j ) {
        if ( j != k ) temp += x(j)*(1+x(j));
      }
      f(k) = x(k)*(2+5*x(k)*x(k)) + 1 - temp;
    }
  }

  integer
  jacobianNnz() const override {
    integer tot = 0;
    for ( integer k = 0; k < n; ++k ) {
      integer k1 = max(0,  k-ml);
      integer k2 = min(n-1,k+mu);
      for ( integer j = k1; j <= k2; ++j ) {
        if ( j != k ) ++tot;
      }
      ++tot;
    }
    return tot;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer k = 0; k < n; ++k ) {
      integer k1 = max(0,  k-ml);
      integer k2 = min(n-1,k+mu);
      for ( integer j = k1; j <= k2; ++j ) {
        if ( j != k ) { ii(kk) = k; jj(kk) = j; ++kk; }
      }
      ii(kk) = jj(kk) = k; ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer k = 0; k < n; ++k ) {
    	integer k1 = max(0,  k-ml);
      integer k2 = min(n-1,k+mu);
      for ( integer j = k1; j <= k2; ++j ) {
        if ( j != k ) jac(kk++) = -(1+2*x(j));
      }
      jac(kk++) = 2+15*x(k)*x(k);
    }
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    for ( integer k = 0; k < n; ++k ) x(k) = -1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
