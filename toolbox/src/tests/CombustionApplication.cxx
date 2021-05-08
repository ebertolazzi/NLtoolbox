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

class CombustionApplication : public nonlinearSystem {

public:

  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  CombustionApplication()
  : nonlinearSystem(
      "Combustion Application",
      "@article{Grosan:2012,\n"
      "  title   = {SOLVING POLYNOMIAL SYSTEMS USING A MODIFIED LINE SEARCH APPROACH},\n"
      "  author  = {Crina Grosan and Ajith Abraham and Vaclav Snasel},\n"
      "  journal = {International Journal of Innovative Computing, Information and Control},\n"
      "  volume  = {8},\n"
      "  number  = {1},\n"
      "  year    = {2012}\n"
      "}\n\n"
      "@book{Morgan:2009,\n"
      "  author = {Morgan, A.},\n"
      "  title  = {Solving Polynomial Systems Using Continuation for\n"
      "            Engineering and Scientific Problems},\n"
      "  publisher = {Society for Industrial and Applied Mathematics},\n"
      "  year = {2009},\n"
      "  doi = {10.1137/1.9780898719031},\n"
      "}\n\n"
      "@article{Hentenryck:1997,\n"
      "  author  = {Van Hentenryck, P. and McAllester, D. and Kapur, D.},\n"
      "  title   = {Solving Polynomial Systems Using a Branch and Prune Approach},\n"
      "  journal = {SIAM Journal on Numerical Analysis},\n"
      "  year    = {1997},\n"
      "  volume  = {34},\n"
      "  number  = {2},\n"
      "  pages   = {797-827},\n"
      "  doi     = {10.1137/S0036142995281504}\n"
      "}\n",
      10
    )
  {
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return x(1) + 2*x(5) + x(8) + 2*x(9) - 1e-5;
      case 1: return x(2) + x(7) - 3e-5;
      case 2: return x(0) + x(2) + 2*x(4) + 2*x(7) + x(8) + x(9) - 5e-5;
      case 3: return x(3) + 2*x(6) - 1e-5;
      case 4: return 0.5140437e-7*x(4) - x(0)*x(0);
      case 5: return 0.1006932e-6*x(5) - 2*x(1)*x(1);
      case 6: return 0.7816278e-15*x(6) - x(3)*x(3);
      case 7: return 0.1496236e-6*x(7) - x(0)*x(2);
      case 8: return 0.6194411e-7*x(8) - x(0)*x(1);
      case 9: return 0.2089296e-14*x(9) - x(0)*x(1);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(1) + 2*x(5) + x(8) + 2*x(9) - 1e-5;
    f(1) = x(2) + x(7) - 3e-5;
    f(2) = x(0) + x(2) + 2*x(4) + 2*x(7) + x(8) + x(9) - 5e-5;
    f(3) = x(3) + 2*x(6) - 1e-5;
    f(4) = 0.5140437e-7*x(4) - x(0)*x(0);
    f(5) = 0.1006932e-6*x(5) - 2*x(1)*x(1);
    f(6) = 0.7816278e-15*x(6) - x(3)*x(3);
    f(7) = 0.1496236e-6*x(7) - x(0)*x(2);
    f(8) = 0.6194411e-7*x(8) - x(0)*x(1);
    f(9) = 0.2089296e-14*x(9) - x(0)*x(1);
  }

  int_type
  jacobianNnz() const override {
    return 29;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk

    SETIJ(0,1); // 1
    SETIJ(0,5);
    SETIJ(0,8);
    SETIJ(0,9);

    SETIJ(1,2); // 5
    SETIJ(1,7);

    SETIJ(2,0); // 7
    SETIJ(2,2);
    SETIJ(2,4);
    SETIJ(2,7);
    SETIJ(2,8);
    SETIJ(2,9);

    SETIJ(3,3); // 13
    SETIJ(3,6);

    SETIJ(4,4); // 14
    SETIJ(4,0);

    SETIJ(5,5); // 17
    SETIJ(5,1);

    SETIJ(6,6); // 19
    SETIJ(6,3);

    SETIJ(7,7); // 21
    SETIJ(7,0);
    SETIJ(7,2);

    SETIJ(8,8); // 24
    SETIJ(8,0);
    SETIJ(8,1);

    SETIJ(9,9); // 27
    SETIJ(9,0);
    SETIJ(9,1);

    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;

    jac(kk++) = 1;
    jac(kk++) = 2;
    jac(kk++) = 1;
    jac(kk++) = 2;

    jac(kk++) = 1;
    jac(kk++) = 1;

    jac(kk++) = 1;
    jac(kk++) = 1;
    jac(kk++) = 2;
    jac(kk++) = 2;
    jac(kk++) = 1;
    jac(kk++) = 1;

    jac(kk++) = 1;
    jac(kk++) = 2;

    jac(kk++) = 0.5140437e-7;
    jac(kk++) = -2*x(0);

    jac(kk++) = 0.1006932e-6;
    jac(kk++) = -4*x(1);

    jac(kk++) = 0.7816278e-15;
    jac(kk++) = -2*x(3);

    jac(kk++) = 0.1496236e-6;
    jac(kk++) = -x(2);
    jac(kk++) = -x(0);

    jac(kk++) = 0.6194411e-7;
    jac(kk++) = -x(1);
    jac(kk++) = -x(0);

    jac(kk++) = 0.2089296e-14;
    jac(kk++) = -x(1);
    jac(kk++) = -x(0);
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(1);
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
