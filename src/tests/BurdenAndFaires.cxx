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

class BurdenAndFaires : public nonlinearSystem {

public:

  BurdenAndFaires()
  : nonlinearSystem(
      "Burden and Faires example 1",
      "@book{burden2005,\n"
      "  author    = {Burden, R. and Faires, J.},\n"
      "  title     = {Numerical Analysis},\n"
      "  year      = {2005},\n"
      "  pages     = {597--640},\n"
      "  publisher = {Thomson Brooks/Cole}\n"
      "}\n",
      3
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return 3*x(0)-cos(x(1)*x(2))-0.5;
      case 1: return x(0)*x(0) - 81 * power2(x(1)+0.1) + sin(x(2)) + 1.06;
      case 2: return exp(-x(0)*x(1)) + 20 * x(2) + (10*m_pi-3)/3;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 3*x(0)-cos(x(1)*x(2))-0.5;
    f(1) = x(0)*x(0) - 81 * power2(x(1)+0.1) + sin(x(2)) + 1.06;
    f(2) = exp(-x(0)*x(1)) + 20 * x(2) + (10*m_pi-3)/3;
  }

  int_type
  jacobianNnz() const override {
    return 9;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = 3;
    jac(1) = sin(x(1)*x(2))*x(2);
    jac(2) = sin(x(1)*x(2))*x(1);

    jac(3) = 2*x(0);
    jac(4) = -162 * (x(1)+0.1);
    jac(5) = cos(x(2));

    jac(6) = -exp(-x(0)*x(1))*x(1);
    jac(7) = -exp(-x(0)*x(1))*x(0);
    jac(8) = 20;
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.1;
    x(1) = 0.1;
    x(2) = -0.1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
