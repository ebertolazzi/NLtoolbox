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

class HelicalValleyFunction : public nonlinearSystem {

public:

  HelicalValleyFunction()
  : nonlinearSystem(
      "Helical valley function",
      "@article{Fletcher:1963,\n"
      "  author  = {Fletcher, R. and Powell, M. J. D.},\n"
      "  title   = {A Rapidly Convergent Descent Method for Minimization},\n"
      "  journal = {The Computer Journal},\n"
      "  year    = {1963},\n"
      "  volume  = {6},\n"
      "  number  = {2},\n"
      "  pages   = {163--168},\n"
      "  doi     = {10.1093/comjnl/6.2.163},\n"
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
      3
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type theta = atan2 ( x(1), x(0) ) / ( 2 * m_pi );
    switch ( k ) {
      case 0: return 10 * ( x(2) - 10 * theta );
      case 1: return 10 * ( sqrt ( power2(x(0)) + power2(x(1)) ) - 1 );
      case 2: return x(2);
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    //real_type theta = atan ( x(1) / x(0) ) / ( 2 * pi );
    //if ( x(0) < 0 ) theta += 0.5;
    real_type theta = atan2 ( x(1), x(0) ) / ( 2 * m_pi );
    f(0) = 10 * ( x(2) - 10 * theta );
    f(1) = 10 * ( sqrt ( power2(x(0)) + power2(x(1)) ) - 1 );
    f(2) = x(2);
  }

  virtual
  int_type
  jacobianNnz() const override
  { return n*n; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type q2 = power2(x(0)) + power2(x(1));
    real_type q  = sqrt ( q2 );
    real_type c  = 50 / m_pi;
    jac(0) =   c * x(1) / q2;
    jac(1) = - c * x(0) / q2;
    jac(2) = 10;

    jac(3) = 10 * x(0) / q;
    jac(4) = 10 * x(1) / q;
    jac(5) = 0;

    jac(6) = 0;
    jac(7) = 0;
    jac(8) = 1;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 0;
    x(2) = 0;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -1;
    x(1) =  0;
    x(2) =  0;
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
