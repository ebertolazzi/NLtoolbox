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

class DennisAndSchnabel2x2example : public nonlinearSystem {
public:
  
  DennisAndSchnabel2x2example()
  : nonlinearSystem(
      "Dennis and Schnabel 2 by 2 example",
      "@book{Dennis:1996,\n"
      "  author    = {Dennis, J. and Schnabel, R.},\n"
      "  title     = {Numerical Methods for Unconstrained\n"
      "               Optimization and Nonlinear Equations},\n"
      "  publisher = {Society for Industrial and Applied Mathematics},\n"
      "  year      = {1996},\n"
      "  doi       = {10.1137/1.9781611971200},\n"
      "}\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return x(0) + x(1) - 3;
      case 1: return power2(x(0)) + power2(x(1))-9;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0) + x(1) - 3;
    f(1) = power2(x(0)) + power2(x(1))-9;
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
    jac(1) = 1;
    jac(2) = 2*x(0);
    jac(3) = 2*x(1);
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 3;
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 5;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
