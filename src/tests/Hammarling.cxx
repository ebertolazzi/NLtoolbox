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

#define HAMMARLING_BIBTEX \
"n by n matrix square root problem (hammarling)" \
"S. J. Hammarling, private communication to P. E. Gill.\n"

class Hammarling2x2matrixSquareRoot : public nonlinearSystem {
  real_type const a00;
  real_type const a01;
  real_type const a10;
  real_type const a11;

public:

  Hammarling2x2matrixSquareRoot()
  : nonlinearSystem(
    "Hammarling 2 by 2 matrix square root problem",
    HAMMARLING_BIBTEX,
    4
  )
  , a00(1.0E-4)
  , a01(1)
  , a10(0)
  , a11(1E-4) {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return (x(0)*x(0) + x(1)*x(2)) - a00;
      case 1: return (x(0)*x(1) + x(1)*x(3)) - a01;
      case 2: return (x(2)*x(0) + x(3)*x(2)) - a10;
      case 3: return (x(2)*x(1) + x(3)*x(3)) - a11;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = (x(0)*x(0) + x(1)*x(2)) - a00;
    f(1) = (x(0)*x(1) + x(1)*x(3)) - a01;
    f(2) = (x(2)*x(0) + x(3)*x(2)) - a10;
    f(3) = (x(2)*x(1) + x(3)*x(3)) - a11;
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
    jac(0)  = 2*x(0);
    jac(1)  = x(2);
    jac(2)  = x(1);
    jac(3)  = 0;
    
    jac(4)  = x(1);
    jac(5)  = x(0) + x(3);
    jac(6)  = 0;
    jac(7)  = x(1);
    
    jac(8)  = x(2);
    jac(9)  = 0;
    jac(10) = x(0) + x(3);
    jac(11) = x(2);
    
    jac(12) = 0;
    jac(13) = x(2);
    jac(14) = x(1);
    jac(15) = 2*x(3);
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1E-2;
    x(1) = 5E+1;
    x(2) = 0;
    x(3) = 1E-2;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 0;
    x(2) = 0;
    x(3) = 1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*
 * Hammarling 3 by 3 matrix square root problem
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblem : public nonlinearSystem {

  real_type a00, a01, a02,
            a10, a11, a12,
            a20, a21, a22;

  dvec_t xe;

  real_type f0( dvec_t const & x ) const
  { return x(0)*x(0) + x(1)*x(3) + x(2)*x(6); }

  real_type f1( dvec_t const & x ) const
  { return x(0)*x(1) + x(1)*x(4) + x(2)*x(7); }

  real_type f2( dvec_t const & x ) const
  { return x(0)*x(2) + x(1)*x(5) + x(2)*x(8); }


  real_type f3( dvec_t const & x ) const
  { return x(3)*x(0) + x(4)*x(3) + x(5)*x(6); }

  real_type f4( dvec_t const & x ) const
  { return x(3)*x(1) + x(4)*x(4) + x(5)*x(7); }

  real_type f5( dvec_t const & x ) const
  { return x(3)*x(2) + x(4)*x(5) + x(5)*x(8); }


  real_type f6( dvec_t const & x ) const
  { return x(6)*x(0) + x(7)*x(3) + x(8)*x(6); }

  real_type f7( dvec_t const & x ) const
  { return x(6)*x(1) + x(7)*x(4) + x(8)*x(7); }

  real_type f8( dvec_t const & x ) const
  { return x(6)*x(2) + x(7)*x(5) + x(8)*x(8); }

public:

  Hammarling3x3matrixSquareRootProblem(
    string const & n,
    real_type x0,
    real_type x1,
    real_type x2,
    real_type x3,
    real_type x4,
    real_type x5,
    real_type x6,
    real_type x7,
    real_type x8
  )
  : nonlinearSystem( n, HAMMARLING_BIBTEX, 9 ) {

    xe.resize(9);

    xe[0] = x0; xe[1] = x1; xe[2] = x2;
    xe[3] = x3; xe[4] = x4; xe[5] = x5;
    xe[6] = x6; xe[7] = x7; xe[8] = x8;
    
    a00 = f0(xe);
    a01 = f1(xe);
    a02 = f2(xe);
    a10 = f3(xe);
    a11 = f4(xe);
    a12 = f5(xe);
    a20 = f6(xe);
    a21 = f7(xe);
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return f0(x) - a00;
      case 1: return f1(x) - a01;
      case 2: return f2(x) - a02;
      case 3: return f3(x) - a10;
      case 4: return f4(x) - a11;
      case 5: return f5(x) - a12;
      case 6: return f6(x) - a20;
      case 7: return f7(x) - a21;
      case 8: return f8(x) - a22;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = f0(x) - a00;
    f(1) = f1(x) - a01;
    f(2) = f2(x) - a02;
    f(3) = f3(x) - a10;
    f(4) = f4(x) - a11;
    f(5) = f5(x) - a12;
    f(6) = f6(x) - a20;
    f(7) = f7(x) - a21;
    f(8) = f8(x) - a22;
  }

  int_type
  jacobianNnz() const override {
    return 45;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0, 0); // 1
    SETIJ(0, 1);
    SETIJ(0, 2);
    SETIJ(0, 3);
    SETIJ(0, 6);

    SETIJ(1, 0); // 6
    SETIJ(1, 1);
    SETIJ(1, 2);
    SETIJ(1, 4);
    SETIJ(1, 7);

    SETIJ(2, 0); // 11
    SETIJ(2, 1);
    SETIJ(2, 2);
    SETIJ(2, 5);
    SETIJ(2, 8);

    SETIJ(3, 0); // 16
    SETIJ(3, 3);
    SETIJ(3, 4);
    SETIJ(3, 5);
    SETIJ(3, 6);

    SETIJ(4, 1); // 21
    SETIJ(4, 3);
    SETIJ(4, 4);
    SETIJ(4, 5);
    SETIJ(4, 7);

    SETIJ(5, 2); // 26
    SETIJ(5, 3);
    SETIJ(5, 4);
    SETIJ(5, 5);
    SETIJ(5, 8);

    SETIJ(6, 0); // 31
    SETIJ(6, 3);
    SETIJ(6, 6);
    SETIJ(6, 7);
    SETIJ(6, 8);

    SETIJ(7, 1); // 36
    SETIJ(7, 4);
    SETIJ(7, 6);
    SETIJ(7, 7);
    SETIJ(7, 8);

    SETIJ(8, 2); // 41
    SETIJ(8, 5);
    SETIJ(8, 6);
    SETIJ(8, 7);
    SETIJ(8, 8);

    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = 2 * x(0);
    jac(kk++) = x(3);
    jac(kk++) = x(6);
    jac(kk++) = x(1);
    jac(kk++) = x(2);

    jac(kk++) = x(1);
    jac(kk++) = x(0) + x(4);
    jac(kk++) = x(7);
    jac(kk++) = x(1);
    jac(kk++) = x(2);

    jac(kk++) = x(2);
    jac(kk++) = x(5);
    jac(kk++) = x(0) + x(8);
    jac(kk++) = x(1);
    jac(kk++) = x(2);
    
    jac(kk++) = x(3);
    jac(kk++) = x(0) + x(4);
    jac(kk++) = x(3);
    jac(kk++) = x(6);
    jac(kk++) = x(5);

    jac(kk++) = x(3);
    jac(kk++) = x(1);
    jac(kk++) = 2*x(4);
    jac(kk++) = x(7);
    jac(kk++) = x(5);

    jac(kk++) = x(3);
    jac(kk++) = x(2);
    jac(kk++) = x(5);
    jac(kk++) = x(4) + x(8);
    jac(kk++) = x(5);
    
    jac(kk++) = x(6);
    jac(kk++) = x(7);
    jac(kk++) = x(0) + x(8);
    jac(kk++) = x(3);
    jac(kk++) = x(6);
    
    jac(kk++) = x(6);
    jac(kk++) = x(7);
    jac(kk++) = x(1);
    jac(kk++) = x(4) + x(8);
    jac(kk++) = x(7);
    
    jac(kk++) = x(6);
    jac(kk++) = x(7);
    jac(kk++) = x(2);
    jac(kk++) = x(5);
    jac(kk++) = 2*x(8);
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x = xe;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 0;
    x(2) = 0;
    x(3) = 0;
    x(4) = 1;
    x(5) = 0;
    x(6) = 0;
    x(7) = 0;
    x(8) = 1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblemN1 : public Hammarling3x3matrixSquareRootProblem {

public:

  Hammarling3x3matrixSquareRootProblemN1()
  : Hammarling3x3matrixSquareRootProblem(
      "Hammarling 3 by 3 matrix square root problem N.1",
      0.01, 50,   0,
      0,    0.01, 0,
      0,    0,    0.01
    )
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblemN2 : public Hammarling3x3matrixSquareRootProblem {

public:

  Hammarling3x3matrixSquareRootProblemN2()
  : Hammarling3x3matrixSquareRootProblem(
      "Hammarling 3 by 3 matrix square root problem N.2",
      0, 0, 1,
      1, 1, 0,
      0, 1, 0
    )
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class Hammarling3x3matrixSquareRootProblemN3 : public Hammarling3x3matrixSquareRootProblem {

public:

  Hammarling3x3matrixSquareRootProblemN3()
  : Hammarling3x3matrixSquareRootProblem(
      "Hammarling 3 by 3 matrix square root problem N.3",
      1, 1, 1,
      0, 0, 0,
      0, 0, 0
    )
  {}

};
