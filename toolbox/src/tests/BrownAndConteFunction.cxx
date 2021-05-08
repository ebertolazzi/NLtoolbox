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

class BrownAndConteFunction : public nonlinearSystem {

  real_type const cst;

public:

  BrownAndConteFunction()
  : nonlinearSystem(
      "Brown and Conte function",
      "@inproceedings{Brown:1967,\n"
      "  author    = {Brown, Kenneth M. and Conte, Samuel D.},\n"
      "  title     = {The Solution of Simultaneous Nonlinear Equations},\n"
      "  booktitle = {Proceedings of the 1967 22Nd National Conference},\n"
      "  series    = {ACM '67},\n"
      "  year      = {1967},\n"
      "  pages     = {111--114},\n"
      "  doi       = {10.1145/800196.805981},\n"
      "  acmid     = {805981},\n"
      "  publisher = {ACM},\n"
      "}\n",
      2
    )
  , cst(1-1/(4*m_pi)) {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
     case 0: return (sin(x(0)*x(1))-x(1)/(2*m_pi)-x(0))/2;
     case 1: return cst*(exp(2*x(0))-m_e)+x(1)*m_e*m_1_pi-2*m_e*x(0);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = (sin(x(0)*x(1))-x(1)/(2*m_pi)-x(0))/2;
    f(1) = cst*(exp(2*x(0))-m_e)+x(1)*m_e*m_1_pi-2*m_e*x(0);
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = (x(1)*cos(x(0)*x(1))-1)/2;
    jac(1) = x(0)*cos(x(0)*x(1))/2-0.25/m_pi;
    jac(2) = 2*(cst*exp(2*x(0))-m_e);
    jac(3) = m_e*m_1_pi;
  }

  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    x(0) = 0.5;
    x(1) = m_pi;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.6;
    x(1) = 3;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
