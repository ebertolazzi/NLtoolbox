/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define HIEBERT_BIBTEX \
"@article{Hiebert:1982,\n" \
"  author  = {Hiebert, K. L.},\n" \
"  title   = {An Evaluation of Mathematical Software That\n" \
"             Solves Systems of Nonlinear Equations},\n" \
"  journal = {ACM Trans. Math. Softw.},\n" \
"  year    = {1982},\n" \
"  volume  = {8},\n" \
"  number  = {1},\n" \
"  pages   = {5--20},\n" \
"  doi     = {10.1145/355984.355986},\n" \
"}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class HiebertChem2x2 : public nonlinearSystem {

public:

  HiebertChem2x2()
  : nonlinearSystem( "Hiebert Chem 2x2", HIEBERT_BIBTEX, 2 )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return x(1) - 10;
      case 1: return x(0)*x(1) -5E4;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(1) - 10;
    f(1) = x(0)*x(1) -5E4;
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
    jac(0) = 0;
    jac(1) = 1;
    jac(2) = x(1);
    jac(3) = x(0);
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 5000;
    x(1) = 10;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 1;
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

class HiebertChem6x6 : public nonlinearSystem {

public:

  HiebertChem6x6()
  : nonlinearSystem( "Hiebert Chem 6x6", HIEBERT_BIBTEX, 6 )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return x(0) + x(1) + x(3) - 0.001;
      case 1: return x(4) + x(5) - 55;
      case 2: return x(0) + x(1) + x(2) + 2 * x(4) + x(5) - 110.001;
      case 3: return x(0) - 0.1 * x(1);
      case 4: return x(0) - 10000 * x(2) * x(3);
      case 5: return x(3) - 55E14 * x(2) * x(5);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0) + x(1) + x(3) - 0.001;
    f(1) = x(4) + x(5) - 55;
    f(2) = x(0) + x(1) + x(2) + 2 * x(4) + x(5) - 110.001;
    f(3) = x(0) - 0.1 * x(1);
    f(4) = x(0) - 10000 * x(2) * x(3);
    f(5) = x(3) - 55E14 * x(2) * x(5);
  }

  int_type
  jacobianNnz() const override {
    return 18;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,3);

    SETIJ(1,4);
    SETIJ(1,5);

    SETIJ(2,0);
    SETIJ(2,1);
    SETIJ(2,2);
    SETIJ(2,4);
    SETIJ(2,5);

    SETIJ(3,0);
    SETIJ(3,1);

    SETIJ(4,0);
    SETIJ(4,2);
    SETIJ(4,3);

    SETIJ(5,2);
    SETIJ(5,3);
    SETIJ(5,5);

    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = 1;
    jac(1) = 1;
    jac(2) = 1;

    jac(3) = 1;
    jac(4) = 1;

    jac(5) = 1;
    jac(6) = 1;
    jac(7) = 1;
    jac(8) = 2;
    jac(9) = 1;

    jac(10) = 1;
    jac(11) = -0.1;

    jac(12) = 1;
    jac(13) = - 10000 * x(3);
    jac(14) = - 10000 * x(2);

    jac(15) = - 55E14 * x(5);
    jac(16) = 1;
    jac(17) = - 55E14 * x(2);
  }

  int_type
  numExactSolution() const override
  { return 3; }

  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
    case 0:
      x << 0.10000000000000865800865802832637712049655126051610e-3,
           0.10000000000000865800865802832637712049655126051610e-2,
           -0.99999999999913419913419799193454109698763804853664e-4,
           -0.10000000000009523809523831159014832546206386567714e-3,
           54.999999999999999818181818181487603305784236699939,
           0.18181818181851239669421576330006082348099514876998e-15;
      break;
    case 1:
      x << -0.18181818182178512396701422689707131359333381750224e-14,
           -0.18181818182178512396701422689707131359333381750224e-13,
           -0.18181818181814876033057918557475582208319103887114e-15,
           0.10000000000200000000003963636363715649586778444953e-2,
           55.001000000000020181818182214512396702144144252600,
          -0.10000000000201818181822145123967021441442526003173e-2;
      break;
    case 2:
      x << 0.82644628099181424635970070847176598631923570438011e-4,
           0.82644628099181424635970070847176598631923570438011e-3,
           0.90909090909186147186147038862875597975247651584838e-4,
           0.90909090909004329004329220681057415048840725181878e-4,
           54.999999999999999818181818182181818181817073593074,
           0.18181818181781818181818292640692640295982173534745e-15;
      break;
    }
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 1;
  }

  int_type
  numInitialPoint() const override
   { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( abs(x(i)) <= 100, "Bad range" );
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class HiebertChem10x10 : public nonlinearSystem {

  real_type const R;
public:

  HiebertChem10x10()
  : nonlinearSystem( "Hiebert Chem 10x10", HIEBERT_BIBTEX, 10)
  , R(10000)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type TOT = x(0)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9);
    switch ( k ) {
      case 0: return x(0)+x(3)-3.0;
      case 1: return 2.0*x(0)+x(1)+x(3)+x(6)+x(7)+x(8)+2.0*x(9)-R;
      case 2: return 2.0*x(1)+2.0*x(4)+x(5)+x(6)-8.0;
      case 3: return 2.0*x(2)+x(4)-4.0*R;
      case 4: return x(0)*x(4)-193.0/1000.0*x(1)*x(3);
      case 5: return x(5)*sqrt(x(1))-2597.0/1000000.0*sqrt(x(1)*x(3)*TOT);
      case 6: return x(6)*sqrt(x(3))-431.0/125000.0*sqrt(x(0)*x(3)*TOT);
      case 7: return x(7)*x(3)-1799.0/100000000.0*x(2)*TOT;
      case 8: return x(8)*x(3)-431.0/2000000.0*x(0)*sqrt(x(2)*TOT);
      case 9: return x(9)*x(3)*x(3)-1923.0/50000000.0*x(3)*x(3)*TOT;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type TOT = x(0)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9);
    f(0) = x(0)+x(3)-3.0;
    f(1) = 2.0*x(0)+x(1)+x(3)+x(6)+x(7)+x(8)+2.0*x(9)-R;
    f(2) = 2.0*x(1)+2.0*x(4)+x(5)+x(6)-8.0;
    f(3) = 2.0*x(2)+x(4)-4.0*R;
    f(4) = x(0)*x(4)-193.0/1000.0*x(1)*x(3);
    f(5) = x(5)*sqrt(x(1))-2597.0/1000000.0*sqrt(x(1)*x(3)*TOT);
    f(6) = x(6)*sqrt(x(3))-431.0/125000.0*sqrt(x(0)*x(3)*TOT);
    f(7) = x(7)*x(3)-1799.0/100000000.0*x(2)*TOT;
    f(8) = x(8)*x(3)-431.0/2000000.0*x(0)*sqrt(x(2)*TOT);
    f(9) = x(9)*x(3)*x(3)-1923.0/50000000.0*x(3)*x(3)*TOT;
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0; // fortran addressing
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type t3 = x(1)*x(3);
    real_type t4 = x(0)+x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9);
    real_type t6 = sqrt(t3*t4);
    real_type t7 = 1/t6;
    real_type t10 = 2597.0/2000000.0*t7*x(1)*x(3);
    real_type t11 = sqrt(x(1));
    real_type t15 = x(3)*t4;
    real_type t25 = x(0)*x(3);
    real_type t27 = sqrt(t25*t4);
    real_type t28 = 1/t27;
    real_type t34 = 431.0/250000.0*t28*x(0)*x(3);
    real_type t35 = sqrt(x(3));
    real_type t45 = 1799.0/100000000.0*x(2);
    real_type t60 = sqrt(x(2)*t4);
    real_type t63 = x(0)/t60;
    real_type t65 = 431.0/4000000.0*t63*x(2);
    real_type t73 = x(3)*x(3);
    real_type t74 = 1923.0/50000000.0*t73;
    jac[caddr(0,0)] = 1.0;
    jac[caddr(0,1)] = 0.0;
    jac[caddr(0,2)] = 0.0;
    jac[caddr(0,3)] = 1.0;
    jac[caddr(0,4)] = 0.0;
    jac[caddr(0,5)] = 0.0;
    jac[caddr(0,6)] = 0.0;
    jac[caddr(0,7)] = 0.0;
    jac[caddr(0,8)] = 0.0;
    jac[caddr(0,9)] = 0.0;
    jac[caddr(1,0)] = 2.0;
    jac[caddr(1,1)] = 1.0;
    jac[caddr(1,2)] = 0.0;
    jac[caddr(1,3)] = 1.0;
    jac[caddr(1,4)] = 0.0;
    jac[caddr(1,5)] = 0.0;
    jac[caddr(1,6)] = 1.0;
    jac[caddr(1,7)] = 1.0;
    jac[caddr(1,8)] = 1.0;
    jac[caddr(1,9)] = 2.0;
    jac[caddr(2,0)] = 0.0;
    jac[caddr(2,1)] = 2.0;
    jac[caddr(2,2)] = 0.0;
    jac[caddr(2,3)] = 0.0;
    jac[caddr(2,4)] = 2.0;
    jac[caddr(2,5)] = 1.0;
    jac[caddr(2,6)] = 1.0;
    jac[caddr(2,7)] = 0.0;
    jac[caddr(2,8)] = 0.0;
    jac[caddr(2,9)] = 0.0;
    jac[caddr(3,0)] = 0.0;
    jac[caddr(3,1)] = 0.0;
    jac[caddr(3,2)] = 2.0;
    jac[caddr(3,3)] = 0.0;
    jac[caddr(3,4)] = 1.0;
    jac[caddr(3,5)] = 0.0;
    jac[caddr(3,6)] = 0.0;
    jac[caddr(3,7)] = 0.0;
    jac[caddr(3,8)] = 0.0;
    jac[caddr(3,9)] = 0.0;
    jac[caddr(4,0)] = x(4);
    jac[caddr(4,1)] = -193.0/1000.0*x(3);
    jac[caddr(4,2)] = 0.0;
    jac[caddr(4,3)] = -193.0/1000.0*x(1);
    jac[caddr(4,4)] = x(0);
    jac[caddr(4,5)] = 0.0;
    jac[caddr(4,6)] = 0.0;
    jac[caddr(4,7)] = 0.0;
    jac[caddr(4,8)] = 0.0;
    jac[caddr(4,9)] = 0.0;
    jac[caddr(5,0)] = -t10;
    jac[caddr(5,1)] = x(5)/t11/2.0-2597.0/2000000.0*t7*(t15+t3);
    jac[caddr(5,2)] = -t10;
    jac[caddr(5,3)] = -2597.0/2000000.0*t7*(x(1)*t4+t3);
    jac[caddr(5,4)] = -t10;
    jac[caddr(5,5)] = t11-t10;
    jac[caddr(5,6)] = -t10;
    jac[caddr(5,7)] = -t10;
    jac[caddr(5,8)] = -t10;
    jac[caddr(5,9)] = -t10;
    jac[caddr(6,0)] = -431.0/250000.0*t28*(t15+t25);
    jac[caddr(6,1)] = -t34;
    jac[caddr(6,2)] = -t34;
    jac[caddr(6,3)] = x(6)/t35/2.0-431.0/250000.0*t28*(x(0)*t4+t25);
    jac[caddr(6,4)] = -t34;
    jac[caddr(6,5)] = -t34;
    jac[caddr(6,6)] = t35-t34;
    jac[caddr(6,7)] = -t34;
    jac[caddr(6,8)] = -t34;
    jac[caddr(6,9)] = -t34;
    jac[caddr(7,0)] = -t45;
    jac[caddr(7,1)] = -t45;
    jac[caddr(7,2)] = -1799.0/100000000.0*x(0)
                      -1799.0/100000000.0*x(1)
                      -1799.0/50000000.0*x(2)
                      -1799.0/100000000.0*x(3)
                      -1799.0/100000000.0*x(4)
                      -1799.0/100000000.0*x(5)
                      -1799.0/100000000.0*x(6)
                      -1799.0/100000000.0*x(7)
                      -1799.0/100000000.0*x(8)
                      -1799.0/100000000.0*x(9);
    jac[caddr(7,3)] = x(7)-t45;
    jac[caddr(7,4)] = -t45;
    jac[caddr(7,5)] = -t45;
    jac[caddr(7,6)] = -t45;
    jac[caddr(7,7)] = x(3)-t45;
    jac[caddr(7,8)] = -t45;
    jac[caddr(7,9)] = -t45;
    jac[caddr(8,0)] = -431.0/2000000.0*t60-t65;
    jac[caddr(8,1)] = -t65;
    jac[caddr(8,2)] = -431.0/4000000.0*t63*(x(0)+x(1)+2.0*x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9));
    jac[caddr(8,3)] = x(8)-t65;
    jac[caddr(8,4)] = -t65;
    jac[caddr(8,5)] = -t65;
    jac[caddr(8,6)] = -t65;
    jac[caddr(8,7)] = -t65;
    jac[caddr(8,8)] = x(3)-t65;
    jac[caddr(8,9)] = -t65;
    jac[caddr(9,0)] = -t74;
    jac[caddr(9,1)] = -t74;
    jac[caddr(9,2)] = -t74;
    jac[caddr(9,3)] = 2.0*x(9)*x(3)-1923.0/25000000.0*t15-t74;
    jac[caddr(9,4)] = -t74;
    jac[caddr(9,5)] = -t74;
    jac[caddr(9,6)] = -t74;
    jac[caddr(9,7)] = -t74;
    jac[caddr(9,8)] = -t74;
    jac[caddr(9,9)] = 49998077.0/50000000.0*t73;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.setZero();
  }

  int_type
  numExactSolution() const override
  { return 3; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill( 1 );
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

