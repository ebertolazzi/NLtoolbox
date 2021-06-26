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

class GriewankFunction : public nonlinearSystem {

public:

  GriewankFunction( integer n )
  : nonlinearSystem(
      "Griewank function",
      "@Article{Griewan:k1981,\n"
      "  author  = {Griewank, A. O.},\n"
      "  title   = {Generalized descent for global optimization},\n"
      "  journal = {Journal of Optimization Theory and Applications},\n"
      "  year    = {1981},\n"
      "  volume  = {34},\n"
      "  number  = {1},\n"
      "  pages   = {11--39},\n"
      "  doi     = {10.1007/BF00933356}\n"
      "}\n",
      n
    )
  {
    UTILS_ASSERT(
      n >= 2 && n <= 20,
      "GriewankFunction(n={} must be in range [2..20]", n
    );
  }

  real_type
  grad( dvec_t const & x, integer k ) const {
    real_type f = 1/sqrt(k+1.0);
    for ( integer i = 0; i < n; ++i ) {
      real_type t = x(i)/sqrt(i+1.0);
      if ( i == k ) f *= sin(t);
      else          f *= cos(t);
    }
    return f;
  }

  real_type
  hess( dvec_t const & x, integer i, integer j ) const {
    if ( i == j ) {
      real_type f = 1/(i+1.0);
      for ( integer k = 0; k < n; ++k ) {
        real_type t = x(k)/sqrt(k+1.0);
        f *= cos(t);
      }
      return f;
    } else {
      real_type f = -1/sqrt((i+1.0)*(j+1.0));
      for ( integer k = 0; k < n; ++k ) {
        real_type t = x(k)/sqrt(k+1.0);
        if ( k == i || k == j ) f *= sin(t);
        else                    f *= cos(t);
      }
      return f;
    }
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    return grad( x, k ) + x(k)/2000.0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const  override{
    for ( integer k = 0; k < n; ++k )
      f(k) = grad( x, k ) + x(k)/2000.0;
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      for ( integer j = 0; j < n; ++j ) {
        jac(kk) = hess(x,i,j);
        if ( i == j ) jac(kk) += 1.0/2000.0;
        ++kk;
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    for ( integer i = 0; i < n; i += 2 ) {
      x(i+0) = -50.0*(1+idx*9);
      x(i+1) =  70.0*(1+idx*9);
    }
  }

  integer
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT0( abs(x(i)) < 1000, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    L.fill(-1000);
    U.fill(1000);
  }

};
