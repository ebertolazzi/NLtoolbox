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

class ExponentialSine : public nonlinearSystem {
public:

  ExponentialSine()
  : nonlinearSystem(
      "Exponential sine",
      "@techreport{Nowak1991,\n"
      "  author = {U. Nowak and L. Weimann},\n"
      "  title  = {A Family of Newton Codes for Systems of Highly Nonlinear Equations},\n"
      "  number = {Technical Report TR-91-10},\n"
      "  year   = {1991}\n"
      "}\n",
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type arg = x(0)*x(0)+x(1)*x(1);
    switch ( k ) {
      case 0: return exp(arg) - 3.0;
      case 1: return x(0)+x(1)-sin(2.0*(x(0)+x(1)));
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type arg = x(0)*x(0)+x(1)*x(1);
    f(0) = exp(arg) - 3.0;
    f(1) = x(0)+x(1)-sin(2.0*(x(0)+x(1)));
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
    real_type arg   = x(0)*x(0)+x(1)*x(1);
    real_type arg_0 = 2*x(0);
    real_type arg_1 = 2*x(1);
    jac(0) = exp(arg)*arg_0;
    jac(1) = exp(arg)*arg_1;
    jac(2) = jac(3) = 1-2*cos(2*(x(0)+x(1)));
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
      case 0: x << 0.81,              0.82;             break;
      case 1: x << 2.99714682530400, -2.95330710241260; break;
      case 2: x << 2.44169014107600, -2.51379895001340; break;
      case 3: x << 2.51913768310200, -2.85296012826840; break;
      case 4: x << -2.26027327648,   -2.41295786268;    break;
    }
  }

  int_type
  numInitialPoint() const override
  { return 5; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
