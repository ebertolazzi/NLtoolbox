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

#define Bohachevsky_BIBTEX \
"@book{Michalewicz:1996,\n" \
"  author = {Michalewicz, Zbigniew},\n" \
"  title = {Genetic Algorithms + Data Structures = Evolution Programs (3rd Ed.)},\n" \
"  year = {1996},\n" \
"  isbn = {3-540-60676-9},\n" \
"  publisher = {Springer-Verlag},\n" \
"  address = {Berlin, Heidelberg},\n" \
"}\n\n" \
"@book{brent2002algorithms,\n" \
"  author={Brent, R.P.},\n" \
"  title={Algorithms for Minimization Without Derivatives},\n" \
"  year={2002},\n" \
"  isbn={9780486419985},\n" \
"  series={Dover Books on Mathematics},\n" \
"  publisher={Dover Publications}\n" \
"}\n"

class BohachevskyN1 : public nonlinearSystem {
public:

  BohachevskyN1()
  : nonlinearSystem("BohachevskyN1",Bohachevsky_BIBTEX,2)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    f(0) = 2.0 * x1 + 0.9 * m_pi * sin( 3.0 * m_pi * x1 );
    f(1) = 4.0 * x2 + 1.6 * m_pi * sin( 4.0 * m_pi * x2 );
  }

  int_type
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    ii(0) = 0; jj(0) = 0;
    ii(1) = 0; jj(1) = 1;
    ii(2) = 1; jj(2) = 0;
    ii(3) = 1; jj(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);

    jac(0) = 2.0 + 2.7 * (m_pi*m_pi) * cos ( 3.0 * m_pi * x1 );
    jac(1) = 0.0;

    jac(2) = 0.0;
    jac(3) = 4.0 + 6.4 * (m_pi*m_pi) * cos ( 4.0 * m_pi * x2 );
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 0;
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
  checkIfAdmissible( dvec_t const & x ) const override {
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BohachevskyN2 : public nonlinearSystem {
public:

  BohachevskyN2()
  : nonlinearSystem("BohachevskyN2",Bohachevsky_BIBTEX,2)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    f(0) = 2 * x1 + 0.9*m_pi*sin(3*m_pi*x1)*cos(4*m_pi*x2);
    f(1) = 4 * x2 + 1.2*m_pi*cos(3*m_pi*x1)*sin(4*m_pi*x2);
  }

  int_type
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    ii(0) = 0; jj(0) = 0;
    ii(1) = 0; jj(1) = 1;
    ii(2) = 1; jj(2) = 0;
    ii(3) = 1; jj(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1  = x(0);
    real_type x2  = x(1);
    real_type pi2 = m_pi*m_pi;
    real_type A   = pi2*sin(3*m_pi*x1)*sin(4*m_pi*x2);
    real_type B   = pi2*cos(3*m_pi*x1)*cos(4*m_pi*x2);
    jac(0) = 2 + 2.7*B;
    jac(1) = jac(2) = -3.6*A;
    jac(3) = 4 + 4.8*B;
  }

  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.6;
    x(1) = 1.3;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BohachevskyN3 : public nonlinearSystem {
public:

  BohachevskyN3()
  : nonlinearSystem("BohachevskyN3",Bohachevsky_BIBTEX,2)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    f(0) = 2.0 * x1 + 0.9 * m_pi * sin( 3*m_pi*x1 );
    f(1) = 4.0 * x2 - 4.0 * m_pi * sin( 4*m_pi*x2 );
  }

  int_type
  jacobianNnz() const override
  { return 2; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    ii(0) = 0; jj(0) = 0;
    ii(1) = 1; jj(1) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type pi2 = m_pi*m_pi;
    jac(0) = 2.0 + 2.7 * pi2*cos(3*m_pi*x1);
    jac(1) = 4.0 - 16.0 * pi2*cos(4*m_pi*x2);
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.5;
    x(1) = 1.0;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
