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

class ShenYpma5 : public nonlinearSystem {
public:

  ShenYpma5()
  : nonlinearSystem(
      "Shen-Ypma Example N.5",
      "@article{,\n"
      "  Author = {Yun-Qiu Shen and Tjalling J. Ypma},\n"
      "  Doi = {10.1016/j.apnum.2004.09.029},\n"
      "  Journal = {Applied Numerical Mathematics},\n"
      "  Number = {2},\n"
      "  Pages = {256 - 265},\n"
      "  Title = {Newton's method for singular nonlinear equations using approximate left and right nullspaces of the Jacobian},\n"
      "  Volume = {54},\n"
      "  Year = {2005},\n"
      "}\n",
      2
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type x0  = x(0);
    real_type x1  = x(1);
    real_type res = 0;
    switch ( i ) {
    case 0: res = x0*x0*(1-x0*x1)+x1*x1; break;
    case 1: res = x0*x0+x1*x1*(3*x0-2);  break;
    }
    return res;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x0 = x(0);
    real_type x1 = x(1);
    f(0) = x0*x0*(1-x0*x1)+x1*x1;
    f(1) = x0*x0+x1*x1*(3*x0-2);
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 4;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(1,0);
    SETIJ(1,1);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type  kk = 0;
    real_type x0 = x(0);
    real_type x1 = x(1);
    jac(kk++) = x0*(2-3*x0*x1);
    jac(kk++) = 2*x1-x0*x0*x0;
    jac(kk++) = 2*x0+3*x1*x1;
    jac(kk++) = x1*(6*x0-4);
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.setZero();
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
    case 0:
      x(0) = x(1) = 0.02;
      break;
    case 1:
      x(0) = 6.2989;
      x(1) = 3.7048;
      break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 2; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ShenYpma7 : public nonlinearSystem {
public:

  ShenYpma7()
  : nonlinearSystem(
      "Shen-Ypma Example N.7",
      "@article{,\n"
      "  Author = {Yun-Qiu Shen and Tjalling J. Ypma},\n"
      "  Doi = {10.1016/j.apnum.2004.09.029},\n"
      "  Journal = {Applied Numerical Mathematics},\n"
      "  Number = {2},\n"
      "  Pages = {256 - 265},\n"
      "  Title = {Newton's method for singular nonlinear equations using approximate left and right nullspaces of the Jacobian},\n"
      "  Volume = {54},\n"
      "  Year = {2005},\n"
      "}\n",
      2
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type x0  = x(0);
    real_type x1  = x(1);
    real_type res = 0;
    switch ( i ) {
    case 0: res = x0*x0-x1*x1;     break;
    case 1: res = 3*(x0*x0-x1*x1); break;
    }
    return res;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x0 = x(0);
    real_type x1 = x(1);
    f(0) = x0*x0-x1*x1;
    f(1) = 3*(x0*x0-x1*x1);
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 4;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(1,0);
    SETIJ(1,1);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type  kk = 0;
    real_type x0 = x(0);
    real_type x1 = x(1);
    jac(kk++) = 2*x0;
    jac(kk++) = -2*x1;
    jac(kk++) = 6*x0;
    jac(kk++) = -6*x1;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.setZero();
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
    case 0:
      x(0) = 0.05;
      x(1) = 0.04;
      break;
    case 1:
      x(0) = 1;
      x(1) = 2;
      break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 2; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ShenYpma8 : public nonlinearSystem {
public:

  ShenYpma8()
  : nonlinearSystem(
      "Shen-Ypma Example N.8",
      "@article{,\n"
      "  Author = {Yun-Qiu Shen and Tjalling J. Ypma},\n"
      "  Doi = {10.1016/j.apnum.2004.09.029},\n"
      "  Journal = {Applied Numerical Mathematics},\n"
      "  Number = {2},\n"
      "  Pages = {256 - 265},\n"
      "  Title = {Newton's method for singular nonlinear equations using approximate left and right nullspaces of the Jacobian},\n"
      "  Volume = {54},\n"
      "  Year = {2005},\n"
      "}\n",
      5
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    real_type x5 = x(4);
    real_type res = 0;
    switch ( i ) {
      case 0: res = x1+x2+x3*x3+x4*x4+x5*x5-2; break;
      case 1: res = x1-x2+x3*x3+x4*x4+x5*x5;    break;
      case 2: res = -x3*x3+x4*x4+x5*x5;         break;
      case 3: res = x3*x3-x4*x4+x5*x5;          break;
      case 4: res = x3*x3+x4*x4-x5*x5;;         break;
    }
    return res;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    real_type x5 = x(4);
    f(0) = x1+x2+x3*x3+x4*x4+x5*x5-2;
    f(1) = x1-x2+x3*x3+x4*x4+x5*x5;
    f(2) = -x3*x3+x4*x4+x5*x5;
    f(3) = x3*x3-x4*x4+x5*x5;
    f(4) = x3*x3+x4*x4-x5*x5;;
  }

  virtual
  int_type
  jacobianNnz() const override {
    return n*n;
  }

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
    //real_type x1 = x(0);
    //real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    real_type x5 = x(4);
    int_type  kk = 0;

    jac(kk++) = 1;
    jac(kk++) = 1;
    jac(kk++) = 2*x3;
    jac(kk++) = 2*x4;
    jac(kk++) = 2*x5;

    jac(kk++) = 1;
    jac(kk++) = -1;
    jac(kk++) = 2*x3;
    jac(kk++) = 2*x4;
    jac(kk++) = 2*x5;

    jac(kk++) = 0;
    jac(kk++) = 0;
    jac(kk++) = -2*x3;
    jac(kk++) = 2*x4;
    jac(kk++) = 2*x5;

    jac(kk++) = 0;
    jac(kk++) = 0;
    jac(kk++) = 2*x3;
    jac(kk++) = -2*x4;
    jac(kk++) = 2*x5;

    jac(kk++) = 0;
    jac(kk++) = 0;
    jac(kk++) = 2*x3;
    jac(kk++) = 2*x4;
    jac(kk++) = -2*x5;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.setZero();
    x(0) = x(1) = 1;
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
    case 0:
      x(0) = 1.02;
      x(1) = 1.02;
      x(2) = 0.02;
      x(3) = 0.02;
      x(4) = 0.02;
      break;
    }
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
