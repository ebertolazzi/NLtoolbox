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

class YixunShi1 : public nonlinearSystem {
public:

  YixunShi1()
  : nonlinearSystem(
      "Shi, Yixun Problem N.3",
      "@article{YixunShi,\n"
      "  Author  = {Shi, Yixun},\n"
      "  Title   = {A globalization procedure for solving nonlinear systems of equations},\n"
      "  Journal = {Numerical Algorithms},\n"
      "  Number  = {2},\n"
      "  Pages   = {273--286},\n"
      "  Volume  = {12},\n"
      "  Year    = {1996},\n"
      "  Doi     = {10.1007/BF02142807},\n"
      "}\n",
      3
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type res = 0;
    switch ( i ) {
    case 0: res = 3*x(0)-cos(x(1)*x(2))-0.5;             break;
    case 1: res = x(0)*x(0)-625*x(1)*x(1);               break;
    case 2: res = exp(-x(0)*x(1))+20*x(2)+(10*m_pi-3)/3; break;
    }
    return res;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 3*x(0)-cos(x(1)*x(2))-0.5;
    f(1) = x(0)*x(0)-625*x(1)*x(1);
    f(2) = exp(-x(0)*x(1))+20*x(2)+(10*m_pi-3)/3;
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 8;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,2);
    SETIJ(1,0);
    SETIJ(1,1);
    SETIJ(2,0);
    SETIJ(2,1);
    SETIJ(2,2);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    real_type x0 = x(0);
    real_type x1 = x(1);
    real_type x2 = x(2);
    real_type t2 = sin(x1*x2);
    real_type t8 = exp(-x0*x1);
    jac(kk++) = 3;
    jac(kk++) = x2*t2;
    jac(kk++) = x1*t2;
    jac(kk++) = 2*x0;
    jac(kk++) = -1250*x1;
    jac(kk++) = -x1*t8;
    jac(kk++) = -x0*t8;
    jac(kk++) = 20;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0.4999816893677071161061094979125234746691;
    x(1) = -0.1999926757470828464424437991650093898676e-1;
    x(2) = -0.5241012469638837054562573412192983462072;
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
    case 0:
      x(0) = x(1) = x(2) = 1;
      break;
    case 1:
      x(0) = 0;
      x(1) = 1e-6;
      x(2) = 0;
      break;
    case 2:
      x(0) = x(1) = x(2) = 0;
      break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 3; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class YixunShi2 : public nonlinearSystem {
public:

  YixunShi2()
  : nonlinearSystem(
      "Shi, Yixun Problem N.4",
      "@article{YixunShi,\n"
      "  Author  = {Shi, Yixun},\n"
      "  Title   = {A globalization procedure for solving nonlinear systems of equations},\n"
      "  Journal = {Numerical Algorithms},\n"
      "  Number  = {2},\n"
      "  Pages   = {273--286},\n"
      "  Volume  = {12},\n"
      "  Year    = {1996},\n"
      "  Doi     = {10.1007/BF02142807},\n"
      "}\n",
      100
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type idx ) const override {
    real_type res = x(idx)*(3-2*x(idx))+1+x(49)/2;
    switch ( idx ) {
    case 0:
      res -= 2*x(1);
      break;
    case 99:
      res -= x(98);
      break;
    default:
      res -= 2*x(idx+1)+x(idx-1);
      break;
    }
    return res;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f = x.array()*(3-2*x.array())+1+x(49)/2;
    f(0)  -= 2*x(1);
    f(99) -= x(98);
    for ( int_type idx = 1; idx < 99; ++idx )
      f(idx) -= 2*x(idx+1)+x(idx-1);
  }

  virtual
  int_type
  jacobianNnz() const override {
    int_type kk = 6;
    for ( int_type idx = 1; idx < 99; ++idx ) {
      kk += 3;
      if ( idx == 48 ) continue;
      if ( idx == 49 ) continue;
      if ( idx == 50 ) continue;
      ++kk;
    }
    return kk;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,49);
    for ( int_type idx = 1; idx < 99; ++idx ) {
      SETIJ(idx,idx-1);
      SETIJ(idx,idx);
      SETIJ(idx,idx+1);
      if ( idx == 48 ) continue;
      if ( idx == 49 ) continue;
      if ( idx == 50 ) continue;
      SETIJ(idx,49);
    }
    SETIJ(99,49);
    SETIJ(99,98);
    SETIJ(99,99);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = 3-4*x(0);
    jac(kk++) = -2;
    jac(kk++) = 0.5;
    for ( int_type idx = 1; idx < 99; ++idx ) {
      jac(kk++) = -1;
      jac(kk++) = 3-4*x(idx);
      jac(kk++) = -2;
      if ( idx == 48 ) {
        jac(kk-1) += 0.5; // (48,49)
      } else if ( idx == 49 ) {
        jac(kk-2) += 0.5; // (49,49)
      } else if ( idx == 50 ) {
        jac(kk-3) += 0.5; // (50,49)
      } else {
        jac(kk++) = 0.5;
      }
    }
    jac(kk++) = 0.5;
    jac(kk++) = -1;
    jac(kk++) = 3-4*x(99);
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    x.resize(n);
    switch ( idx ) {
    case 0:
      x.fill(-1);
      break;
    case 1:
      x.fill(0.04);
      break;
    case 2:
      x.fill(0.044285);
      break;
    case 3:
      x.fill(1);
      break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 4; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class YixunShi3 : public nonlinearSystem {
public:

  YixunShi3()
  : nonlinearSystem(
      "Shi, Yixun Problem N.5",
      "@article{YixunShi,\n"
      "  Author  = {Shi, Yixun},\n"
      "  Title   = {A globalization procedure for solving nonlinear systems of equations},\n"
      "  Journal = {Numerical Algorithms},\n"
      "  Number  = {2},\n"
      "  Pages   = {273--286},\n"
      "  Volume  = {12},\n"
      "  Year    = {1996},\n"
      "  Doi     = {10.1007/BF02142807},\n"
      "}\n",
      100
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type idx ) const override {
    real_type res = 1-x(99)+0.5*x(98)-x(97)-x(96)+3*x(95)+x(idx)*(3-2*x(idx));
    switch ( idx ) {
    case 0:
      res -= 2*x(1);
      break;
    case 99:
      res -= x(98);
      break;
    default:
      res -= 2*x(idx+1)+x(idx-1);
      break;
    }
    return res;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f = 1-x(99)+0.5*x(98)-x(97)-x(96)+3*x(95)+x.array()*(3-2*x.array());
    f(0)  -= 2*x(1);
    f(99) -= x(98);
    for ( int_type idx = 1; idx < 99; ++idx )
      f(idx) -= 2*x(idx+1)+x(idx-1);
  }

  virtual
  int_type
  jacobianNnz() const override {
    int_type kk = 0;
    for ( int_type idx = 0; idx < 100; ++idx ) {
      kk += 5;
      if ( idx == 0 ) {
        kk += 2;
      } else if ( idx == 94 ) {
        kk += 2;
      } else if ( idx == 95 ) {
        ++kk;
      } else if ( idx < 94 ) {
        kk += 3;
      }
    }
    return kk;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type idx = 0; idx < 100; ++idx ) {
      SETIJ(idx,95);
      SETIJ(idx,96);
      SETIJ(idx,97);
      SETIJ(idx,98);
      SETIJ(idx,99);
      if ( idx == 0 ) {
        SETIJ(0,0);
        SETIJ(0,1);
      } else if ( idx == 94 ) {
        SETIJ(94,93);
        SETIJ(94,94);
      } else if ( idx == 95 ) {
        SETIJ(95,94);
      } else if ( idx < 94 ) {
        SETIJ(idx,idx-1);
        SETIJ(idx,idx);
        SETIJ(idx,idx+1);
      }
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type idx = 0; idx < 100; ++idx ) {
      jac(kk++) = 3;    // 95
      jac(kk++) = -1;   // 96
      jac(kk++) = -1;   // 97
      jac(kk++) = 0.5;  // 98
      jac(kk++) = -1;   // 99
      if ( idx == 0 ) {
        jac(kk++) = 3-4*x(0);
        jac(kk++) = -2;
      } else if ( idx < 94 ) {
        jac(kk++) = -1;
        jac(kk++) = 3-4*x(idx);
        jac(kk++) = -2;
      } else if ( idx == 94 ) {
        jac(kk-5) -= 2; // (94,95)
        jac(kk++) = -1;
        jac(kk++) = 3-4*x(idx);
      } else if ( idx == 95 ) {
        jac(kk-5) += 3-4*x(idx); // (95,95)
        jac(kk-4) -= 2;// (95,96)
        jac(kk++) = -1;
      } else if ( idx == 96 ) {
        jac(kk-5) -= 1;
        jac(kk-4) += 3-4*x(idx);
        jac(kk-3) -= 2;
      } else if ( idx == 97 ) {
        jac(kk-4) -= 1;
        jac(kk-3) += 3-4*x(idx);
        jac(kk-2) -= 2;
      } else if ( idx == 98 ) {
        jac(kk-3) -= 1;
        jac(kk-2) += 3-4*x(idx);
        jac(kk-1) -= 2;
      } else if ( idx == 99 ) {
        jac(kk-2) -= 1;
        jac(kk-1) += 3-4*x(idx);
      }
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    x.resize(n);
    switch ( idx ) {
    case 0:
      x.fill(-1);
      break;
    case 1:
      x.fill(0.04);
      break;
    case 2:
      x.fill(0.043283);
      break;
    case 3:
      x.fill(1);
      break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 4; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class YixunShi4 : public nonlinearSystem {
  real_type h;
public:

  YixunShi4()
  : nonlinearSystem(
      "Shi, Yixun Problem N.6 (Singular Broyden)",
      "@article{YixunShi,\n"
      "  Author  = {Shi, Yixun},\n"
      "  Title   = {A globalization procedure for solving nonlinear systems of equations},\n"
      "  Journal = {Numerical Algorithms},\n"
      "  Number  = {2},\n"
      "  Pages   = {273--286},\n"
      "  Volume  = {12},\n"
      "  Year    = {1996},\n"
      "  Doi     = {10.1007/BF02142807},\n"
      "}\n",
      100
    )
  , h(0.5)
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type idx ) const override {
    real_type res = (3-h*x(idx))*x(idx)+1;
    switch ( idx ) {
    case 0:
      res -= 2*x(1);
      break;
    case 99:
      res -= x(98);
      break;
    default:
      res -= 2*x(idx+1)+x(idx-1);
      break;
    }
    return res*res;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f      = (3-h*x.array())*x.array()+1;
    f(0)  -= 2*x(1);
    f(99) -= x(98);
    for ( int_type idx = 1; idx < 99; ++idx )
      f(idx) -= 2*x(idx+1)+x(idx-1);
    f = f.array()*f.array();
  }

  virtual
  int_type
  jacobianNnz() const override {
    int_type kk = 4;
    for ( int_type idx = 1; idx < 99; ++idx ) kk += 3;
    return kk;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    for ( int_type idx = 1; idx < 99; ++idx ) {
      SETIJ(idx,idx-1);
      SETIJ(idx,idx);
      SETIJ(idx,idx+1);
    }
    SETIJ(99,98);
    SETIJ(99,99);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    real_type tmp   = (3-h*x(0))*x(0)-2*x(1)+1;
    real_type tmp_1 = 3-2*h*x(0);
    jac(kk++) = 2*tmp*tmp_1;
    jac(kk++) = -4*tmp;
    for ( int_type idx = 1; idx < 99; ++idx ) {
      tmp   = (3-h*x(idx))*x(idx)-x(idx-1)-2*x(idx+1)+1;
      tmp_1 = 3-2*h*x(idx);
      jac(kk++) = -2*tmp;
      jac(kk++) = 2*tmp*tmp_1;
      jac(kk++) = -4*tmp;
    }
    tmp   = (3-h*x(99))*x(99)-x(98)+1;
    tmp_1 = 3-2*h*x(99);
    jac(kk++) = -2*tmp;
    jac(kk++) = 2*tmp*tmp_1;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    x.resize(n);
    switch ( idx ) {
    case 0:
      x.fill(-1);
      break;
    case 1:
      x.fill(-0.001);
      break;
    case 2:
      x.setZero();
      break;
    case 3:
      x.fill(0.1);
      break;
    case 4:
      x.fill(0.17);
      break;
    case 5:
      x.fill(1);
      break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 6; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
