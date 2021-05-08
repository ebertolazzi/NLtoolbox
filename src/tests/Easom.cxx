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

class Easom : public nonlinearSystem {
public:

  Easom()
  : nonlinearSystem(
      "Easom Function",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      2
    )
  {
  }

  real_type
  evalFk ( dvec_t const & x, int_type k ) const override {
    real_type arg = - power2(x(0)-m_pi) - power2(x(1)-m_pi);
    switch ( k ) {
    case 0: return cos(x(1))*(sin(x(0)) + 2*cos(x(0))*(x(0)-m_pi))*exp(arg);
    case 1: return cos(x(0))*(sin(x(1)) + 2*cos(x(1))*(x(1)-m_pi))*exp(arg);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type arg = - power2(x(0)-m_pi) - power2(x(1)-m_pi);
    f(0) = cos(x(1))*(sin(x(0)) + 2*cos(x(0))*(x(0)-m_pi))*exp(arg);
    f(1) = cos(x(0))*(sin(x(1)) + 2*cos(x(1))*(x(1)-m_pi))*exp(arg);
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
    real_type arg     = - power2(x(0)-m_pi) - power2(x(1)-m_pi);
    real_type dargdx1 = -2 * (x(0) - m_pi );
    real_type dargdx2 = -2 * (x(1) - m_pi );

    real_type factor = cos(x(1)) * (sin(x(0)) - cos(x(0)) * dargdx1 );
    real_type dfdx1  = cos(x(1)) * (cos(x(0)) + sin(x(0)) * dargdx1 + 2 * cos(x(0)) );
    real_type dfdx2  = -sin(x(1)) * (sin(x(0)) - cos(x(0)) * dargdx1 );

    jac(0) = ( dfdx1 + factor * dargdx1 ) * exp ( arg );
    jac(1) = ( dfdx2 + factor * dargdx2 ) * exp ( arg );

    factor = cos(x(0)) * (sin(x(1)) - cos(x(1)) * dargdx2 );
    dfdx1  = -sin(x(0)) * (sin(x(1)) - cos(x(1)) * dargdx2 );
    dfdx2  = cos(x(0)) * (cos(x(1)) + sin(x(1)) * dargdx2 + 2*cos(x(1)) );

    jac(2) = ( dfdx1 + factor * dargdx1 ) * exp ( arg );
    jac(3) = ( dfdx2 + factor * dargdx2 ) * exp ( arg );
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = m_pi;
    x(1) = m_pi;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.5;
    x(1) = 1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
