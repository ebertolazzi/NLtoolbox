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

class BoxProblem : public nonlinearSystem {
  real_type const tmp;

public:

  BoxProblem()
  : nonlinearSystem(
      "Box Problem",
      "@article{Box:1966,\n"
      "  author  = {Box, M. J.},\n"
      "  title   = {A Comparison of Several Current Optimization Methods,\n"
      "             and the use of Transformations in Constrained Problems},\n"
      "  journal = {The Computer Journal},\n"
      "  volume  = {9},\n"
      "  number  = {1},\n"
      "  pages   = {67-77},\n"
      "  year    = {1966},\n"
      "  doi     = {10.1093/comjnl/9.1.67},\n"
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
      3
    )
  , tmp( exp(-0.1) + exp(-1) )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    switch ( k ) {
      case 0: return exp(-0.1*x(0))-exp(-0.1*x(1))-x(2)*tmp;
      case 1: return exp(-0.2*x(0))-exp(-0.2*x(1))-x(2)*tmp;
      case 2: return exp(-0.3*x(0))-exp(-0.3*x(1))-x(2)*tmp;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = exp(-0.1*x(0))-exp(-0.1*x(1))-x(2)*tmp;
    f(1) = exp(-0.2*x(0))-exp(-0.2*x(1))-x(2)*tmp;
    f(2) = exp(-0.3*x(0))-exp(-0.3*x(1))-x(2)*tmp;
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
    jac(0) = -0.1*exp(-0.1*x(0));
    jac(1) =  0.1*exp(-0.1*x(1));
    jac(2) = -tmp;

    jac(3) = -0.2*exp(-0.2*x(0));
    jac(4) =  0.2*exp(-0.2*x(1));
    jac(5) = -tmp;

    jac(6) = -0.3*exp(-0.3*x(0));
    jac(7) =  0.3*exp(-0.3*x(1));
    jac(8) = -tmp;
  }

  integer
  numExactSolution() const override
  { return 2; }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
    case 0:
      x(0) = 1;
      x(1) = 10;
      x(2) = 1;
      break;
    case 1:
      x(0) =  1;
      x(1) =  10;
      x(2) = -1;
      break;
    }
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 1;
    x(2) = 0;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
