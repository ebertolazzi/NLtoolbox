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

class ChemicalEquilibriumApplication : public nonlinearSystem {
  real_type const R;
  real_type const R5;
  real_type const R6;
  real_type const R7;
  real_type const R8;
  real_type const R9;
  real_type const R10;
public:

  ChemicalEquilibriumApplication()
  : nonlinearSystem(
      "Chemical Equilibrium Application",
      "@article{Meintjes:1990,\n"
      "  author  = {Meintjes, Keith and Morgan, Alexander P.},\n"
      "  title   = {Chemical Equilibrium Systems As Numerical Test Problems},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1990},\n"
      "  volume  = {16},\n"
      "  number  = {2},\n"
      "  pages   = {143--151},\n"
      "  doi     = {10.1145/78928.78930},\n"
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
      5
    )
  , R(10)
  , R5(0.193)
  , R6(0.002597/sqrt(40.0))
  , R7(0.003448/sqrt(40.0))
  , R8(0.00001799/sqrt(40.0))
  , R9(0.0002155/sqrt(40.0))
  , R10(0.00003846/sqrt(40.0))
  {
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
    case 0:
      return x(0)*x(1) + x(0) - 3*x(4);
    case 1:
      return 2*x(0)*x(1) + x(0) + x(1)*x(2)*x(2)
             + R8*x(1) - R*x(4) + 2*R10*x(1)*x(1)
             + R7*x(1)*x(2) + R9*x(1)*x(3);
    case 2:
      return 2*(x(1)+R5)*x(2)*x(2) - 8*x(4) + R6*x(2) + R7*x(1)*x(2);
    case 3:
      return (R9*x(1)+2*x(3))*x(3)-4*R*x(4);
    case 4:
      return x(0)*(x(1)+1) + (R8+R10*x(1)+R7*x(2)+R9*x(3))*x(1) +
             (R5+x(1))*x(2)*x(2) + x(3)*x(3)-1 + R6*x(2);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0)*x(1) + x(0) - 3*x(4);
    f(1) = 2*x(0)*x(1) + x(0) + x(1)*x(2)*x(2)
         + R8*x(1)-R*x(4) + 2*R10*x(1)*x(1)
         + R7*x(1)*x(2) + R9*x(1)*x(3);
    f(2) = 2*(x(1)+R5)*x(2)*x(2) - 8*x(4) + R6*x(2) + R7*x(1)*x(2);
    f(3) = (R9*x(1)+2*x(3))*x(3) - 4*R*x(4);
    f(4) = x(0)*(x(1)+1) + (R8+R10*x(1)+R7*x(2)+R9*x(3))*x(1)
         + (R5+x(1))*x(2)*x(2) + x(3)*x(3)-1 + R6*x(2);
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
    SETIJ(0,4); // 3

    SETIJ(1,0);
    SETIJ(1,1);
    SETIJ(1,2);
    SETIJ(1,3);
    SETIJ(1,4); // 5

    SETIJ(2,1);
    SETIJ(2,2);
    SETIJ(2,4); // 3

    SETIJ(3,1);
    SETIJ(3,3);
    SETIJ(3,4); // 3

    SETIJ(4,0);
    SETIJ(4,1);
    SETIJ(4,2);
    SETIJ(4,3); // 4

    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;

    jac(kk++) = 1+x(1);
    jac(kk++) = x(0);
    jac(kk++) = -3;

    jac(kk++) = 2*x(1) + 1;
    jac(kk++) = 4*R10*x(1)+R7*x(2)+R9*x(3)+x(2)*x(2)+R8+2*x(0);
    jac(kk++) = 2*x(1)*x(2) + R7*x(1);
    jac(kk++) = R9*x(1);
    jac(kk++) = -R;

    jac(kk++) = 2*x(2)*x(2)+R7*x(2);
    jac(kk++) = 4*(x(1)+R5)*x(2)+R6+R7*x(1);
    jac(kk++) = -8;

    jac(kk++) = R9*x(3);
    jac(kk++) = R9*x(1)+ 4*x(3);
    jac(kk++) = -4*R;

    jac(kk++) = x(1)+1;
    jac(kk++) = x(0) + R8+2*R10*x(1)+R7*x(2)+R9*x(3) + x(2)*x(2);
    jac(kk++) = R7*x(1) + 2*(R5+x(1))*x(2) + R6;
    jac(kk++) = R9*x(1) + 2*x(3);

  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override { return 0; }

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
