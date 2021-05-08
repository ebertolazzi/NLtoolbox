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

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    return x(0) * power2(x(0)-5);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0) * power2(x(0)-5);
  }

  virtual
  int_type
  jacobianNnz() const override
  { return 1; }

  virtual
  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = (3*x(0)-5)*(x(0)-5);
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0; break;
      case 1: x(0) = 5; break;
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 2; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
