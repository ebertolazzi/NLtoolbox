/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

// TEST 222

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BroydenTridiagonalFunction : public nonlinearSystem {
  real_type const alpha;
  real_type const beta;
public:
  
  BroydenTridiagonalFunction( real_type _alpha,
                              real_type _beta,
                              int_type  _neq )
  : nonlinearSystem(
      "Broyden tridiagonal function",
      "@article{Broyden:1965,\n"
      "  author  = {Broyden, C. G.},\n"
      "  title   = {A class of methods for solving nonlinear simultaneous equations},\n"
      "  journal = {Mathematics of Computation},\n"
      "  volume  = {19},\n"
      "  year    = {1965},\n"
      "  pages   = {577--593},\n"
      "  doi     = {10.2307/2003941}\n"
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
      _neq
    )
  , alpha(_alpha)
  , beta(_beta)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = (3-alpha*x(k))*x(k)+beta;
    if ( k > 0   ) f -= x(k-1);
    if ( k < n-1 ) f -= 2*x(k+1);
    return f;
  }
  
  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k ) {
    	f(k) = (3-alpha*x(k))*x(k)+beta;
      if ( k > 0   ) f(k) -= x(k-1);
      if ( k < n-1 ) f(k) -= 2*x(k+1);
    }
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
    for ( int_type k = 0; k < n;   ++k ) { SETIJ(k,k); }
    for ( int_type k = 1; k < n;   ++k ) { SETIJ(k,k-1); }
    for ( int_type k = 0; k < n-1; ++k ) { SETIJ(k,k+1); }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n;   ++k ) jac(kk++) = 3 - 2*alpha*x(k);
    for ( int_type k = 0; k < n-1; ++k ) jac(kk++) = -1;
    for ( int_type k = 0; k < n-1; ++k ) jac(kk++) = -2;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k ) x(k) = -1;
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
