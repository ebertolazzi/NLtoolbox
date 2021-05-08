
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

#define BUNLSI_BIBTEX \
"@article{Buzzi:1986,\n" \
"  author  = {Guido Buzzi Ferraris and Enrico Tronconi},\n" \
"  title   = {Bunlsi—{A} fortran program for solution of systems\n" \
"             of nonlinear algebraic equations},\n" \
"  journal = {Computers \\& Chemical Engineering},\n" \
"  volume  = {10},\n" \
"  number  = {2},\n" \
"  pages   = {129--141},\n" \
"  year    = {1986},\n" \
"  doi     = {10.1016/0098-1354(86)85025-6},\n" \
"}\n"

class BUNLSI5 : public nonlinearSystem {

public:

  BUNLSI5()
  : nonlinearSystem( "BUNLSI example 5", BUNLSI_BIBTEX, 8 )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    real_type x5 = x(4);
    real_type x6 = x(5);
    real_type x7 = x(6);
    real_type x8 = x(7);
    switch ( k ) {
     case 0: return x1-1;
     case 1: return x2-sqrt(x1) - exp(x1) - 15;
     case 2: return x3-x1*x2/100 - sin(x1) - 1;
     case 3: return x4 - power2(x2+x1+x3/2) + 150;
     case 4: return x5 - (x4-x2)/(x1*x2*x3);
     case 5: return x6 - x5*pow(x1,1.0/3.0) -exp(x5);
     case 6: return x7 - (x1 - sqrt(x4) - x5*x5)*exp(x5);
     case 7: return x8 - x1 - x5 - x6 - x3;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    real_type x5 = x(4);
    real_type x6 = x(5);
    real_type x7 = x(6);
    real_type x8 = x(7);
    f(0) = x1-1;
    f(1) = x2 - sqrt(x1) - exp(x1) - 15;
    f(2) = x3 - x1*x2/100 - sin(x1) - 1;
    f(3) = x4 - power2(x2+x1+x3/2) + 150;
    f(4) = x5 - (x4-x2)/(x1*x2*x3);
    f(5) = x6 - x5*pow(x1,1.0/3.0) -exp(x5);
    f(6) = x7 - (x1 - sqrt(x4) - x5*x5)*exp(x5);
    f(7) = x8 - x1 - x5 - x6 - x3;
  }

  int_type
  jacobianNnz() const override {
    return 27;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I-1; jj(kk) = J-1; ++kk
    SETIJ(1,1); // 1

    SETIJ(2,1);
    SETIJ(2,2);

    SETIJ(3,1); // 4
    SETIJ(3,2);
    SETIJ(3,3);

    SETIJ(4,1); // 7
    SETIJ(4,2);
    SETIJ(4,3);
    SETIJ(4,4);

    SETIJ(5,1); // 11
    SETIJ(5,2);
    SETIJ(5,3);
    SETIJ(5,4);
    SETIJ(5,5);

    SETIJ(6,1); // 16
    SETIJ(6,5);
    SETIJ(6,6);

    SETIJ(7,1); // 19
    SETIJ(7,4);
    SETIJ(7,5);
    SETIJ(7,7);

    SETIJ(8,1); // 23
    SETIJ(8,3);
    SETIJ(8,5);
    SETIJ(8,6);
    SETIJ(8,8);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    real_type x5 = x(4);
    //real_type x6 = x(5);
    //real_type x7 = x(6);
    //real_type x8 = x(7);

    int_type kk = 0;
    jac(kk++) = 1;

    jac(kk++) = -0.5/sqrt(x1)-exp(x1);
    jac(kk++) = 1;

    jac(kk++) = -x2/100-cos(x1);
    jac(kk++) = -x1/100;
    jac(kk++) = 1;

    jac(kk++) = -2*(x2+x1+x3/2);
    jac(kk++) = -2*(x2+x1+x3/2);
    jac(kk++) = -(x2+x1+x3/2);
    jac(kk++) = 1;

    jac(kk++) = (x4-x2)/(x1*x1*x2*x3);
    jac(kk++) = x4/(x1*x2*x2*x3);
    jac(kk++) = (x4-x2)/(x1*x2*x3*x3);
    jac(kk++) = -1/(x1*x2*x3);
    jac(kk++) = 1;

    jac(kk++) = -x5*pow(x1,-2.0/3.0)/3;
    jac(kk++) = -pow(x1,1.0/3.0)-exp(x5);
    jac(kk++) = 1;

    jac(kk++) = -exp(x5);
    jac(kk++) = 0.5*exp(x5)/sqrt(x4);
    jac(kk++) = (x5*(x5+2)+sqrt(x4)-x1)*exp(x5);
    jac(kk++) = 1;

    jac(kk++) = -1;
    jac(kk++) = -1;
    jac(kk++) = -1;
    jac(kk++) = -1;
    jac(kk++) =  1;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 18.71;
    x(2) = 2.0286;
    x(3) = 279.8;
    x(4) = 6.876;
    x(5) = 976.15;
    x(6) = -61079.6;
    x(7) = 986.06;
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 20;
    x(2) = 2.2;
    x(3) = 100;
    x(4) = 2;
    x(5) = 8;
    x(6) = -60;
    x(7) = 15;
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

class BUNLSI6 : public nonlinearSystem {

public:

  BUNLSI6()
  : nonlinearSystem( "BUNLSI example 6", BUNLSI_BIBTEX, 30 )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type x1  = x(0);
    real_type x2  = x(1);
    real_type x3  = x(2);
    real_type x4  = x(3);
    real_type x5  = x(4);
    real_type x6  = x(5);
    real_type x7  = x(6);
    real_type x8  = x(7);
    real_type x9  = x(8);
    real_type x10 = x(9);
    real_type x11 = x(10);
    real_type x12 = x(11);
    real_type x13 = x(12);
    real_type x14 = x(13);
    real_type x15 = x(14);
    real_type x16 = x(15);
    real_type x17 = x(16);
    real_type x18 = x(17);
    real_type x19 = x(18);
    real_type x20 = x(19);
    real_type x21 = x(20);
    real_type x22 = x(21);
    real_type x23 = x(22);
    real_type x24 = x(23);
    real_type x25 = x(24);
    real_type x26 = x(25);
    real_type x27 = x(26);
    real_type x28 = x(27);
    real_type x29 = x(28);
    real_type x30 = x(29);
    switch ( k ) {
      case 0:  return x1 - 1;
      case 1:  return 6*x1 - x2;
      case 2:  return 5.4*x1 + x2 - x3;
      case 3:  return x2 + x3 - x4;
      case 4:  return 100*(x1+x2-x5);
      case 5:  return x5 + 100*x4 - x6;
      case 6:  return x1 - 0.1*x4 - x7;
      case 7:  return x5 + 90*x2 - x8;
      case 8:  return x1 + x7 - x9;
      case 9:  return x9 + x3 - x10;
      case 10: return x1*x2 - x11;
      case 11: return x11/x3 - x12;
      case 12: return sqrt(x5) + x4 - x13;
      case 13: return log(x1*x6) + x10 - x14;
      case 14: return -sin(x2) + log(x2+x10) - x15;
      case 15: return x3*x9*x10 - x14 - x16;
      case 16: return 69.1*x1 - 0.01*x2*x5 - x17;
      case 17: return x6/(x5*x15) - x18;
      case 18: return sqrt(x10*x16) - x19;
      case 19: return pow(x1,2.6) + x12*x12 - x20;
      case 20: return x11*exp(x20) - x21;
      case 21: return x5/x10 + x13 - x22;
      case 22: return 0.1*x2 + x1*x10/x4 - x23;
      case 23: return log(x5*x8) + x2 - x24;
      case 24: return x11*x24 - x22 - x25;
      case 25: return power2(x11*x12*x14) - x26 - 1500;
      case 26: return x1*x5/10000 - x27;
      case 27: return power2(x6-x8)/x22 - x28;
      case 28: return x16*x19/x21 - x29;
      case 29: return x27*x28 + x22 + sqrt(x28) - x30;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1  = x(0);
    real_type x2  = x(1);
    real_type x3  = x(2);
    real_type x4  = x(3);
    real_type x5  = x(4);
    real_type x6  = x(5);
    real_type x7  = x(6);
    real_type x8  = x(7);
    real_type x9  = x(8);
    real_type x10 = x(9);
    real_type x11 = x(10);
    real_type x12 = x(11);
    real_type x13 = x(12);
    real_type x14 = x(13);
    real_type x15 = x(14);
    real_type x16 = x(15);
    real_type x17 = x(16);
    real_type x18 = x(17);
    real_type x19 = x(18);
    real_type x20 = x(19);
    real_type x21 = x(20);
    real_type x22 = x(21);
    real_type x23 = x(22);
    real_type x24 = x(23);
    real_type x25 = x(24);
    real_type x26 = x(25);
    real_type x27 = x(26);
    real_type x28 = x(27);
    real_type x29 = x(28);
    real_type x30 = x(29);
    f(0)  = x1 - 1;
    f(1)  = 6*x1 - x2;
    f(2)  = 5.4*x1 + x2 - x3;
    f(3)  = x2 + x3 - x4;
    f(4)  = 100*(x1+x2-x5);
    f(5)  = x5 + 100*x4 - x6;
    f(6)  = x1 - 0.1*x4 - x7;
    f(7)  = x5 + 90*x2 - x8;
    f(8)  = x1 + x7 - x9;
    f(9)  = x9 + x3 - x10;
    f(10) = x1*x2 - x11;
    f(11) = x11/x3 - x12;
    f(12) = sqrt(x5) + x4 - x13;
    f(13) = log(x1*x6) + x10 - x14;
    f(14) = -sin(x2) + log(x2+x10) - x15;
    f(15) = x3*x9*x10 - x14 - x16;
    f(16) = 69.1*x1 - 0.01*x2*x5 - x17;
    f(17) = x6/(x5*x15) - x18;
    f(18) = sqrt(x10*x16) - x19;
    f(19) = pow(x1,2.6) + x12*x12 - x20;
    f(20) = x11*exp(x20) - x21;
    f(21) = x5/x10 + x13 - x22;
    f(22) = 0.1*x2 + x1*x10/x4 - x23;
    f(23) = log(x5*x8) + x2 - x24;
    f(24) = x11*x24 - x22 - x25;
    f(25) = power2(x11*x12*x14) - x26 - 1500;
    f(26) = x1*x5/10000 - x27;
    f(27) = power2(x6-x8)/x22 - x28;
    f(28) = x16*x19/x21 - x29;
    f(29) = x27*x28 + x22 + sqrt(x28) - x30;
  }

  int_type
  jacobianNnz() const override {
    return 101;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I-1; jj(kk) = J-1; ++kk

    SETIJ(1,1); // +1

    SETIJ(2,1); // 1
    SETIJ(2,2);

    SETIJ(3,1); // 3
    SETIJ(3,2);
    SETIJ(3,3);

    SETIJ(4,2); // 6
    SETIJ(4,3);
    SETIJ(4,4);

    SETIJ(5,1); // 9
    SETIJ(5,2);
    SETIJ(5,5);

    SETIJ(6,4); // 12
    SETIJ(6,5);
    SETIJ(6,6);

    SETIJ(7,1); // 15
    SETIJ(7,4);
    SETIJ(7,7);

    SETIJ(8,2); // 18
    SETIJ(8,5);
    SETIJ(8,8);

    SETIJ(9,1); // 21
    SETIJ(9,7);
    SETIJ(9,9);

    SETIJ(10,3); // 24
    SETIJ(10,9);
    SETIJ(10,10);

    SETIJ(11,1); // 27
    SETIJ(11,2);
    SETIJ(11,11);

    SETIJ(12,3); // 30
    SETIJ(12,11);
    SETIJ(12,12);

    SETIJ(13,4); // 33
    SETIJ(13,5);
    SETIJ(13,13);

    SETIJ(14,1); // 36
    SETIJ(14,6);
    SETIJ(14,10);
    SETIJ(14,14);

    SETIJ(15,2); // 40
    SETIJ(15,10);
    SETIJ(15,15);

    SETIJ(16,3); // 43
    SETIJ(16,9);
    SETIJ(16,10);
    SETIJ(16,14);
    SETIJ(16,16);

    SETIJ(17,1); // 48
    SETIJ(17,2);
    SETIJ(17,5);
    SETIJ(17,17);

    SETIJ(18,5); // 52
    SETIJ(18,6);
    SETIJ(18,15);
    SETIJ(18,18);

    SETIJ(19,10); // 56
    SETIJ(19,16);
    SETIJ(19,19);

    SETIJ(20,1); // 59
    SETIJ(20,12);
    SETIJ(20,20);

    SETIJ(21,11); // 62
    SETIJ(21,20);
    SETIJ(21,21);

    SETIJ(22,5); // 65
    SETIJ(22,10);
    SETIJ(22,13);
    SETIJ(22,22);

    SETIJ(23,1); // 69
    SETIJ(23,2);
    SETIJ(23,4);
    SETIJ(23,10);
    SETIJ(23,23);

    SETIJ(24,2); // 74
    SETIJ(24,5);
    SETIJ(24,8);
    SETIJ(24,24);

    SETIJ(25,11); // 78
    SETIJ(25,22);
    SETIJ(25,24);
    SETIJ(25,25);

    SETIJ(26,11); // 82
    SETIJ(26,12);
    SETIJ(26,14);
    SETIJ(26,26);

    SETIJ(27,1); // 86
    SETIJ(27,5);
    SETIJ(27,27);

    SETIJ(28,6); // 89
    SETIJ(28,8);
    SETIJ(28,22);
    SETIJ(28,28);

    SETIJ(29,16); // 93
    SETIJ(29,19);
    SETIJ(29,21);
    SETIJ(29,29);

    SETIJ(30,22); // 97
    SETIJ(30,27);
    SETIJ(30,28);
    SETIJ(30,30);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    real_type x5 = x(4);
    real_type x6 = x(5);
    //real_type x7 = x(6);
    real_type x8 = x(7);
    real_type x9 = x(8);
    real_type x10 = x(9);
    real_type x11 = x(10);
    real_type x12 = x(11);
    //real_type x13 = x(12);
    real_type x14 = x(13);
    real_type x15 = x(14);
    real_type x16 = x(15);
    //real_type x17 = x(16);
    //real_type x18 = x(17);
    real_type x19 = x(18);
    real_type x20 = x(19);
    real_type x21 = x(20);
    real_type x22 = x(21);
    //real_type x23 = x(22);
    real_type x24 = x(23);
    //real_type x25 = x(24);
    //real_type x26 = x(25);
    real_type x27 = x(26);
    real_type x28 = x(27);
    //real_type x29 = x(28);
    //real_type x30 = x(29);

    int_type kk = 0;

    jac(kk++) = 1;

    jac(kk++) = 6;
    jac(kk++) = -1;

    jac(kk++) = 5.4;
    jac(kk++) = 1;
    jac(kk++) = -1;

    jac(kk++) = 1;
    jac(kk++) = 1;
    jac(kk++) = -1;

    jac(kk++) = 100;
    jac(kk++) = 100;
    jac(kk++) = -100;

    jac(kk++) = 100;
    jac(kk++) = 1;
    jac(kk++) = -1;

    jac(kk++) = 1;
    jac(kk++) = -0.1;
    jac(kk++) = -1;

    jac(kk++) = 90;
    jac(kk++) = 1;
    jac(kk++) = -1;

    jac(kk++) = 1;
    jac(kk++) = 1;
    jac(kk++) = -1;

    jac(kk++) = 1;
    jac(kk++) = 1;
    jac(kk++) = -1;

    jac(kk++) = x2;
    jac(kk++) = x1;
    jac(kk++) = -1;

    jac(kk++) = -x11/(x3*x3);
    jac(kk++) = 1/x3;
    jac(kk++) = -1;

    jac(kk++) = 1;
    jac(kk++) = 0.5/sqrt(x5);
    jac(kk++) = -1;

    jac(kk++) = 1/x1;
    jac(kk++) = 1/x6;
    jac(kk++) = 1;
    jac(kk++) = -1;

    jac(kk++) = 1/(x2+x10)-cos(x2);
    jac(kk++) = 1/(x2+x10);
    jac(kk++) = -1;

    jac(kk++) = x9*x10;
    jac(kk++) = x3*x10;
    jac(kk++) = x3*x9;
    jac(kk++) = -1;
    jac(kk++) = -1;

    jac(kk++) = 69.1;
    jac(kk++) = -0.01*x5;
    jac(kk++) = -0.01*x2;
    jac(kk++) = -1;

    jac(kk++) = -x6/(x5*x5*x15);
    jac(kk++) = 1/(x5*x15);
    jac(kk++) = -x6/(x5*x15*x15);
    jac(kk++) = -1;

    jac(kk++) = 0.5*x16/sqrt(x10*x16);
    jac(kk++) = 0.5*x10/sqrt(x10*x16);
    jac(kk++) = -1;

    jac(kk++) = 2.6*pow(x1,1.6);
    jac(kk++) = 2*x12;
    jac(kk++) = -1;

    jac(kk++) = exp(x20);
    jac(kk++) = x11*exp(x20);
    jac(kk++) = -1;

    jac(kk++) = 1/x10;
    jac(kk++) = -x5/(x10*x10);
    jac(kk++) =  1;
    jac(kk++) = -1;

    jac(kk++) = x10/x4;
    jac(kk++) = 0.1;
    jac(kk++) = -x1*x10/(x4*x4);
    jac(kk++) = x1/x4;
    jac(kk++) = -1;

    jac(kk++) = 1;
    jac(kk++) = 1/x5;
    jac(kk++) = 1/x8;
    jac(kk++) = -1;

    jac(kk++) = x24;
    jac(kk++) = -1;
    jac(kk++) = x11;
    jac(kk++) = -1;

    jac(kk++) = 2*(x11*x12*x14)*x12*x14;
    jac(kk++) = 2*(x11*x12*x14)*x11*x14;
    jac(kk++) = 2*(x11*x12*x14)*x11*x12;
    jac(kk++) = -1;

    jac(kk++) = x5/10000;
    jac(kk++) = x1/10000;
    jac(kk++) = -1;

    jac(kk++) = 2*(x6-x8)/x22;
    jac(kk++) = -2*(x6-x8)/x22;
    jac(kk++) = -power2((x6-x8)/x22);
    jac(kk++) = -1;

    jac(kk++) = x19/x21;
    jac(kk++) = x16/x21;
    jac(kk++) = -x16*x19/power2(x21);
    jac(kk++) = -1;

    jac(kk++) = 1;
    jac(kk++) = x28;
    jac(kk++) = x27+0.5/sqrt(x28);
    jac(kk++) = -1;
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0)  = 1.5;
    x(1)  = 1;
    x(2)  = 50;
    x(3)  = 50;
    x(4)  = 500;
    x(5)  = 5000;
    x(6)  = 0;
    x(7)  = 5000;
    x(8)  = 0;
    x(9)  = 5000;
    x(10) = 0.5;
    x(11) = 50;
    x(12) = 5;
    x(13) = 0.5;
    x(14) = 50;
    x(15) = 50;
    x(16) = 5;
    x(17) = 50;
    x(18) = 50;
    x(19) = 5;
    x(20) = 50;
    x(21) = 500;
    x(22) = 5;
    x(23) = 50;
    x(24) = -50;
    x(25) = 5000;
    x(26) = 0.05;
    x(27) = 50000;
    x(28) = 5;
    x(29) = 5000;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    real_type x1  = x(0);
    real_type x4  = x(3);
    real_type x5  = x(4);
    real_type x6  = x(5);
    real_type x8  = x(7);
    real_type x10 = x(9);
    real_type x22 = x(21);
    real_type x28 = x(27);
    NONLIN_ASSERT( x1  >  0, "x1" );
    NONLIN_ASSERT( x4  >  0, "x4" );
    NONLIN_ASSERT( x5  >  0, "x5" );
    NONLIN_ASSERT( x6  >  0, "x6" );
    NONLIN_ASSERT( x8  >  0, "x8" );
    NONLIN_ASSERT( x10 >  0, "x10" );
    NONLIN_ASSERT( x22 >  0, "x22" );
    NONLIN_ASSERT( x28 >  0, "x28" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(real_max);
    L.fill(-real_max);
    L[0] = L[3] = L[4] = L[5] = L[7] = L[9] = L[21] = L[27] = 0;
  }
};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/
