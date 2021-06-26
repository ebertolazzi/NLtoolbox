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

class ExtendedEigerSikorskiStenger : public nonlinearSystem {
public:

  ExtendedEigerSikorskiStenger()
  : nonlinearSystem(
      "Extended Eiger-Sikorski-Stenger Function",
      "@article{Eiger:1984,\n"
      "  author  = {Eiger, A. and Sikorski, K. and Stenger, F.},\n"
      "  title   = {A Bisection Method for Systems of Nonlinear Equations},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1984},\n"
      "  volume  = {10},\n"
      "  number  = {4},\n"
      "  pages   = {367--377},\n"
      "  doi     = {10.1145/2701.2705},\n"
      "}\n\n"
      "@article{Kearfott:1987,\n"
      "  author  = {Kearfott, R. Baker},\n"
      "  title   = {Some Tests of Generalized Bisection},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1987},\n"
      "  volume  = {13},\n"
      "  number  = {3},\n"
      "  pages   = {197--220},\n"
      "  doi     = {10.1145/29380.29862},\n"
      "}\n",
      9
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( i == n-1 ) return power2(x(n-1)-0.1) + x(0) - 0.1;
    return power2( x(i) - 0.1) + x(i+1) - 0.1;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n-1; ++i )
      f(i) = power2( x(i) - 0.1) + x(i+1) - 0.1;
    f(n-1) = power2(x(n-1)-0.1) + x(0) - 0.1;
  }

  integer
  jacobianNnz() const override
  { return 2*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n-1; ++i ) {
      ii(kk) = jj(kk) = i; ++kk;
      ii(kk) = i; jj(kk) = i+1; ++kk;
    }
    ii(kk) = jj(kk) = n-1; ++kk;
    ii(kk) = n-1; jj(kk) = 0;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n-1; ++i ) {
      jac(kk++) = 2*(x(i) - 0.1);
      jac(kk++) = 1;
    }
    jac(kk++) = 2*(x(n-1)-0.1);
    jac(kk++) = 1;
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.fill(-2000);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
