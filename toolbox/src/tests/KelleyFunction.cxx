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

class KelleyFunction : public nonlinearSystem {

public:

  KelleyFunction()
  : nonlinearSystem(
      "Kelley Function",
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
  evalFk( dvec_t const & x_in, integer k ) const override {
    real_type x = x_in[0];
    real_type y = x_in[1];
    if ( k == 0 ) return x*x+y*y-2;
    return exp(x-1)+y*y-2;
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x = x_in[0];
    real_type y = x_in[1];
    f(0) = x*x+y*y-2;
    f(1) = exp(x-1)+y*y-2;
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
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x = x_in[0];
    real_type y = x_in[1];
    jac(0) = 2*x;
    jac(1) = 2*y;
    jac(2) = exp(x-1);
    jac(3) = 2*y;
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 1;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 2;
    x(1) = 1e-6;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
