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

class DarvishiBarati : public nonlinearSystem {
public:

  DarvishiBarati()
  : nonlinearSystem(
      "DarvishiBarati",
      "@article{Darvishi:2007,\n"
      "  author  = {Darvishi, M.T. and Barati, A.},\n"
      "  title   = {Super cubic iterative methods to solve systems\n"
      "             of nonlinear equations},\n"
      "  journal = {Applied Mathematics and Computation},\n"
      "  volume  = {188},\n"
      "  number  = {2},\n"
      "  pages   = {1678--1685},\n"
      "  year    = {2007},\n"
      "  doi     = {10.1016/j.amc.2006.11.022}\n"
      "}\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    switch ( k ) {
    case 0: return exp(x1+x2)+x1*cos(x2);
    case 1: return x1+x2-1;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    f(0) = exp(x1+x2)+x1*cos(x2);
    f(1) = x1+x2-1;
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    jac(0) = exp(x1+x2)+cos(x2);
    jac(1) = exp(x1+x2)-x1*sin(x2);
    jac(2) = 1;
    jac(3) = 1;
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = -4;
    x(1) = 5;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
