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

class McCormicFunction : public nonlinearSystem {
public:

  McCormicFunction()
  : nonlinearSystem(
      "McCormic function",
      "MVF - Multivariate Test Functions Library in C for Unconstrained Global Optimization\n"
      "Ernesto P. Adorio Department of Mathematics U.P. Diliman\n"
      "ernesto.adorio@gmail.com eadorio@yahoo.com, 2005\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type t2 = cos(x(0)+x(1));
    real_type t3 = 2.0*x(0);
    real_type t4 = 2.0*x(1);
    switch ( k ) {
      case 0: return t2+t3-t4-3.0/2.0;
      case 1: return t2-t3+t4+5.0/2.0;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type t2 = cos(x(0)+x(1));
    real_type t3 = 2.0*x(0);
    real_type t4 = 2.0*x(1);
    f(0) = t2+t3-t4-3.0/2.0;
    f(1) = t2-t3+t4+5.0/2.0;
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
    real_type t2 = sin(x(0)+x(1));
    real_type t3 = -t2+2.0;
    real_type t4 = -t2-2.0;
    jac(0) = t3;
    jac(1) = t4;
    jac(2) = t4;
    jac(3) = t3;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = m_pi/3+0.5;
    x(1) = m_pi/3-0.5;
  }

  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
      case 0 : x(0) = -1.5; x(1) = 4; break;
      case 1 : x(0) = -3;   x(1) = 4; break;
    }
  }

  int_type
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
