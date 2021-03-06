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

class SampleProblem19 : public nonlinearSystem {
public:
  
  SampleProblem19()
  : nonlinearSystem(
      "Sample problem 19",
      "no doc",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    switch ( k ) {
      case 0: return x(0) * ( power2(x(0)) + power2(x(1)) );
      case 1: return x(1) * ( power2(x(0)) + power2(x(1)) );
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0) * ( power2(x(0)) + power2(x(1)) );
    f(1) = x(1) * ( power2(x(0)) + power2(x(1)) );
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
    jac(0) = 3*power2(x(0)) + power2(x(1));
    jac(1) = 2*x(0)*x(1);

    jac(2) = 2*x(0)*x(1);
    jac(3) = power2(x(0))+3*power2(x(1));
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 3;
    x(1) = 3;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
