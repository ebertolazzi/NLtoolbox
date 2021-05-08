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

class FreudensteinRothFunction : public nonlinearSystem {
public:

  FreudensteinRothFunction()
  : nonlinearSystem(
      "Freudenstein-Roth function",
      "@article{Freudenstein:1963,\n"
      "  author  = {Freudenstein, Ferdinand and Roth, Bernhard},\n"
      "  title   = {Numerical Solution of Systems of Nonlinear Equations},\n"
      "  journal = {J. ACM},\n"
      "  year    = {1963},\n"
      "  volume  = {10},\n"
      "  number  = {4},\n"
      "  pages   = {550--556},\n"
      "  doi     = {10.1145/321186.321200}\n"
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
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return x(0)-power3(x(1))+5*power2(x(1))-2*x(1)-13;
      case 1: return x(0)+power3(x(1))+power2(x(1))-14*x(1)-29;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0)-power3(x(1))+5*power2(x(1))-2*x(1)-13;
    f(1) = x(0)+power3(x(1))+power2(x(1))-14*x(1)-29;
  }

  int_type
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
    jac(0) = 1;
    jac(1) = -3*power2(x(1))+10*x(1)-2;
    jac(2) = 1;
    jac(3) = 3*power2(x(1))+2*x(1)-14;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 5;
    x(1) = 4;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) =  0.5;
    x(1) = -2;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
