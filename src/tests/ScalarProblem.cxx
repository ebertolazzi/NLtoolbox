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

class ScalarProblem : public nonlinearSystem {
public:

  ScalarProblem()
  : nonlinearSystem(
      "Scalar problem f(x) = x * ( x - 5 )**2",
      "no doc",
      1
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    return x(0) * power2(x(0)-5);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0) * power2(x(0)-5);
  }

  integer
  jacobianNnz() const override
  { return 1; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = (3*x(0)-5)*(x(0)-5);
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0; break;
      case 1: x(0) = 5; break;
    }
  }

  integer
  numExactSolution() const override
  { return 2; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
