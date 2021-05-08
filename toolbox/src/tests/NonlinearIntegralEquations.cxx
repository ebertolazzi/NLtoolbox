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

class NonlinearIntegralEquations : public nonlinearSystem {

public:

  NonlinearIntegralEquations()
  : nonlinearSystem(
      "Nonlinear Integral Equations",
      "@article{More:1979,\n"
      "  author  = {Mor{\'e}, Jorge J. and Cosnard, Michel Y.},\n"
      "  title   = {Numerical Solution of Nonlinear Equations},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1979},\n"
      "  volume  = {5},\n"
      "  number  = {1},\n"
      "  pages   = {64--85},\n"
      "  doi     = {10.1145/355815.355820},\n"
      "}\n",
      100
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type  k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type h = 1.0/(n-1.0);
    
    for ( int_type i = 0; i < n; ++i ) f(i) = x(i);

    real_type si = 0;
    for ( int_type i = 1; i < n-1; ++i ) {
      real_type t = h*i;
      si   += power3(x(i)+t+1);
      f(i) += 0.5*(1-t)*si;
    }

    si = 0;
    for ( int_type i = n-2; i > 0; --i ) {
      real_type t = h*i;
      f(i) += 0.5*(1-t)*t*si;
      si   += power2(x(i)+t+1);
    }
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0; // fortran addressing
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {

    jac.setZero();

    real_type h = 1.0/(n-1.0);

    for ( int_type i = 0; i < n; ++i ) jac[caddr(i,i)] = 1;

    for ( int_type i = 0; i < n; ++i ) {
      real_type t = h*i;
      for ( int_type j = 1; j <= i; ++j ) {
        jac[caddr(i,j)] += 1.5*(1-t)*power2(x(j)+h*j+1);
      }
      for ( int_type j = i+1; j < n-1; ++j ) {
        jac[caddr(i,j)] += (1-t)*t*(x(j)+h*j+1);
      }
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
