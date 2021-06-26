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

class BraninRCOS : public nonlinearSystem {
  real_type const a, d, e, b, c, ff;
public:

  BraninRCOS()
  : nonlinearSystem(
      "BraninRCOS",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      2
    )
  , a(1.0)
  , d(6.0)
  , e(10.0)
  , b(5.1/(4.0 * m_pi* m_pi))
  , c(5.0/m_pi)
  , ff(1.0/(8.0*m_pi))
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    f(0) = 2.0*a*(x2 - b*x1*x1 + c * x1 - d ) * (c-2*b*x1) - e*(1-ff)*sin(x1);
    f(1) = 2.0*a*(x2 - b*x1*x1 + c * x1 - d );
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
    real_type x1 = x(0);
    real_type x2 = x(1);

    jac(0) = 2*a*power2(c-2*b*x1)
           - 4*a*b*(x2-b*x1*x1+c*x1-d)
           - e*(1-ff)*cos(x1);
    jac(1) = jac(2) = 2*a*(c-2*b*x1);
    jac(3) = 2*a;
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) =   -m_pi; x(1) = 12.275; break;
      case 1: x(0) =    m_pi; x(1) =  2.275; break;
      case 2: x(0) = 9.42478; x(1) =  2.475; break;
    }
  }

  integer
  numExactSolution() const override
  { return 3; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = -1;
    x(1) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
