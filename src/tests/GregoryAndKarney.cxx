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

class GregoryAndKarney : public nonlinearSystem {
public:

  GregoryAndKarney( integer n )
  : nonlinearSystem(
      "Gregory and Karney Tridiagonal Matrix Function",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      n
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    if ( k == 0 ) {
      return x(0)-x(1)-2;
    } else if ( k == n-1 ) {
      return 2*x(k)-x(k-1);
    } else {
      return 2*x(k)-x(k-1)-x(k+1);
    }
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n; ++i ) {
      if ( i == 0 ) {
        f(0) = x(0)-x(1)-2;
      } else if ( i==n-1 ) {
        f(i) = 2*x(i)-x(i-1);
      } else {
        f(i) = 2*x(i)-x(i-1)-x(i+1);
      }
    }
  }

  integer
  jacobianNnz() const override {
    return 3*n-2;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( integer i = 0; i < n; ++i ) {
      if ( i == 0 ) {
        SETIJ(i,i);
        SETIJ(i,i+1);
      } else if ( i==n-1 ) {
        SETIJ(i,i);
        SETIJ(i,i-1);
      } else {
        SETIJ(i,i-1);
        SETIJ(i,i);
        SETIJ(i,i+1);
      }
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      if ( i == 0 ) {
        jac(kk++) = 1;
        jac(kk++) = -1;
      } else if ( i==n-1 ) {
        jac(kk++) = 2;
        jac(kk++) = -1;
      } else {
        jac(kk++) = -1;
        jac(kk++) = 2;
        jac(kk++) = -1;
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    for ( integer i = 0; i < n; ++i ) x(i) = n-i;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.setZero();
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  UTILS_ASSERT( abs(x(i)) < 10000, "x[" << i << "] = "<< x(i) << " too big");
  }

};
