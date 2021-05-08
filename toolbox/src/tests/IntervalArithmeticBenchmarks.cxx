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

class IntervalArithmeticBenchmarks : public nonlinearSystem {
  real_type a1, a2, a3, a4, a5, a6, a7, a8, a9, a10;
  real_type b1, b2, b3, b4, b5, b6, b7, b8, b9, b10;
public:

  IntervalArithmeticBenchmarks()
  : nonlinearSystem(
      "Interval Arithmetic Benchmarks",
      "@article{Morgan:1987,\n"
      "  author  = {Alexander Morgan and Andrew Sommese},\n"
      "  title   = {Computing all solutions to polynomial\n"
      "             systems using homotopy continuation},\n"
      "  journal = {Applied Mathematics and Computation},\n"
      "  volume  = {24},\n"
      "  number  = {2},\n"
      "  pages   = {115--138},\n"
      "  year    = {1987},\n"
      "  issn    = {0096-3003},\n"
      "  doi     = {10.1016/0096-3003(87)90064-6}\n"
      "}\n\n"
      "@article{Hentenryck:1997,\n"
      "  author  = {Van Hentenryck, P. and McAllester, D. and Kapur, D.},\n"
      "  title   = {Solving Polynomial Systems Using a Branch and Prune Approach},\n"
      "  journal = {SIAM Journal on Numerical Analysis},\n"
      "  year    = {1997},\n"
      "  volume  = {34},\n"
      "  number  = {2},\n"
      "  pages   = {797-827},\n"
      "  doi = {10.1137/S0036142995281504}\n"
      "}\n",
      10
    )
  {
    a1  = 0.25428722;
    a2  = 0.37842197;
    a3  = 0.27162577;
    a4  = 0.19807914;
    a5  = 0.44166728;
    a6  = 0.14654113;
    a7  = 0.42937161;
    a8  = 0.07056438;
    a9  = 0.34504906;
    a10 = 0.42651102;

    b1  = 0.18324757;
    b2  = 0.16275449;
    b3  = 0.16955071;
    b4  = 0.15585316;
    b5  = 0.19950920;
    b6  = 0.18922793;
    b7  = 0.21180486;
    b8  = 0.17081208;
    b9  = 0.19612740;
    b10 = 0.21466544;
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return x(0) - a1  - b1  * x(3)*x(2)*x(8);
      case 1: return x(1) - a2  - b2  * x(0)*x(9)*x(5);
      case 2: return x(2) - a3  - b3  * x(0)*x(1)*x(9);
      case 3: return x(3) - a4  - b4  * x(6)*x(0)*x(5);
      case 4: return x(4) - a5  - b5  * x(6)*x(5)*x(2);
      case 5: return x(5) - a6  - b6  * x(7)*x(4)*x(9);
      case 6: return x(6) - a7  - b7  * x(1)*x(4)*x(7);
      case 7: return x(7) - a8  - b8  * x(0)*x(6)*x(5);
      case 8: return x(8) - a9  - b9  * x(9)*x(5)*x(7);
      case 9: return x(9) - a10 - b10 * x(3)*x(7)*x(0);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0) - a1  - b1  * x(3)*x(2)*x(8);
    f(1) = x(1) - a2  - b2  * x(0)*x(9)*x(5);
    f(2) = x(2) - a3  - b3  * x(0)*x(1)*x(9);
    f(3) = x(3) - a4  - b4  * x(6)*x(0)*x(5);
    f(4) = x(4) - a5  - b5  * x(6)*x(5)*x(2);
    f(5) = x(5) - a6  - b6  * x(7)*x(4)*x(9);
    f(6) = x(6) - a7  - b7  * x(1)*x(4)*x(7);
    f(7) = x(7) - a8  - b8  * x(0)*x(6)*x(5);
    f(8) = x(8) - a9  - b9  * x(9)*x(5)*x(7);
    f(9) = x(9) - a10 - b10 * x(3)*x(7)*x(0);
  }

  int_type
  jacobianNnz() const override {
    return 40;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk

    SETIJ(0,0);
    SETIJ(0,2);
    SETIJ(0,3);
    SETIJ(0,8);

    SETIJ(1,1);
    SETIJ(1,0);
    SETIJ(1,5);
    SETIJ(1,9);

    SETIJ(2,2);
    SETIJ(2,0);
    SETIJ(2,1);
    SETIJ(2,9);

    SETIJ(3,3);
    SETIJ(3,0);
    SETIJ(3,5);
    SETIJ(3,6);

    SETIJ(4,4);
    SETIJ(4,2);
    SETIJ(4,5);
    SETIJ(4,6);

    SETIJ(5,5);
    SETIJ(5,4);
    SETIJ(5,7);
    SETIJ(5,9);

    SETIJ(6,6);
    SETIJ(6,1);
    SETIJ(6,4);
    SETIJ(6,7);

    SETIJ(7,7);
    SETIJ(7,0);
    SETIJ(7,6);
    SETIJ(7,5);

    SETIJ(8,8);
    SETIJ(8,5);
    SETIJ(8,7);
    SETIJ(8,9);

    SETIJ(9,9);
    SETIJ(9,3);
    SETIJ(9,7);
    SETIJ(9,0);

    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;

    jac(kk++) = 1;
    jac(kk++) = -b1*x(3)*x(8);
    jac(kk++) = -b1*x(2)*x(8);
    jac(kk++) = -b1*x(3)*x(2);

    jac(kk++) = 1;
    jac(kk++) = -b2*x(9)*x(5);
    jac(kk++) = -b2*x(0)*x(9);
    jac(kk++) = -b2*x(0)*x(5);

    jac(kk++) = 1;
    jac(kk++) = -b3*x(1)*x(9);
    jac(kk++) = -b3*x(0)*x(9);
    jac(kk++) = -b3*x(0)*x(1);

    jac(kk++) = 1;
    jac(kk++) = -b4*x(6)*x(5);
    jac(kk++) = -b4*x(6)*x(0);
    jac(kk++) = -b4*x(0)*x(5);

    jac(kk++) = 1;
    jac(kk++) = -b5*x(6)*x(5);
    jac(kk++) = -b5*x(6)*x(2);
    jac(kk++) = -b5*x(5)*x(2);

    jac(kk++) = 1;
    jac(kk++) = -b6*x(7)*x(9);
    jac(kk++) = -b6*x(4)*x(9);
    jac(kk++) = -b6*x(7)*x(4);

    jac(kk++) = 1;
    jac(kk++) = -b7*x(4)*x(7);
    jac(kk++) = -b7*x(1)*x(7);
    jac(kk++) = -b7*x(1)*x(4);

    jac(kk++) = 1;
    jac(kk++) = -b8*x(6)*x(5);
    jac(kk++) = -b8*x(0)*x(5);
    jac(kk++) = -b8*x(0)*x(6);

    jac(kk++) = 1;
    jac(kk++) = -b9*x(9)*x(7);
    jac(kk++) = -b9*x(9)*x(5);
    jac(kk++) = -b9*x(5)*x(7);

    jac(kk++) = 1;
    jac(kk++) = -b10*x(7)*x(0);
    jac(kk++) = -b10*x(3)*x(0);
    jac(kk++) = -b10*x(3)*x(7);

  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.setZero();
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
