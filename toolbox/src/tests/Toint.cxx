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

class Toint225 : public nonlinearSystem {
public:
  
  Toint225( int_type neq )
  : nonlinearSystem(
      "Toint N.225",
      "@Article{Spedicato1997,\n"
      "  author  = {Spedicato, E. and Huang, Z.},\n"
      "  title   = {Numerical experience with newton-like methods\n"
      "             for nonlinear algebraic systems},\n"
      "  journal = {Computing},\n"
      "  year    = {1997},\n"
      "  volume  = {58},\n"
      "  number  = {1},\n"
      "  pages   = {69--89},\n"
      "  doi     = {10.1007/BF02684472},\n"
      "}\n",
      neq
    )
  { checkMinEquations(n,2); }

  real_type
  phi1( real_type s, real_type t ) const
  { return 3*s*s+2*t-5+sin(s-t)*sin(s+t); }

  real_type
  phi2( real_type s, real_type t ) const
  { return 4*t-3+s*exp(s-t); }

  real_type
  phi1_1( real_type s, real_type t ) const
  { return 6*s + cos(s - t)*sin(s + t) + sin(s - t)*cos(s + t); }
  
  real_type
  phi1_2( real_type s, real_type t ) const
  { return 2 - cos(s - t)*sin(s + t) + sin(s - t)*cos(s + t); }

  real_type
  phi2_1( real_type s, real_type t ) const
  { return (1+s)*exp(s - t); }

  real_type
  phi2_2( real_type s, real_type t ) const
  { return 4 - s*exp(s - t); }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    if      ( k == 0   ) f = phi1(x(0),x(1));
    else if ( k == n-1 ) f = phi2(x(n-2),x(n-1));
    else                 f = phi1(x(k),x(k+1))+phi2(x(k-1),x(k));
    return f;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = phi1(x(0),x(1));
    for ( int_type k = 1; k < n-1; ++k )
      f(k) = phi1(x(k),x(k+1))+phi2(x(k-1),x(k));
    f(n-1) = phi2(x(n-2),x(n-1));
  }

  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    for ( int_type k = 1; k < n-1; ++k ) {
      SETIJ(k,k-1);
      SETIJ(k,k);
      SETIJ(k,k+1);
    }
    SETIJ(n-1,n-2);
    SETIJ(n-1,n-1);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = phi1_1(x(0),x(1));
    jac(kk++) = phi1_2(x(0),x(1));
    for ( int_type k = 1; k < n-1; ++k ) {
      jac(kk++) = phi2_1(x(k-1),x(k));
      jac(kk++) = phi1_1(x(k),x(k+1))+phi2_2(x(k-1),x(k));
      jac(kk++) = phi1_2(x(k),x(k+1));
    }
    jac(kk++) = phi2_1(x(n-2),x(n-1));
    jac(kk++) = phi2_2(x(n-2),x(n-1));
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }
  
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k ) x(k) = 0;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
