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

class CubeFunction : public nonlinearSystem {
public:

  CubeFunction()
  : nonlinearSystem(
      "Cube Function",
      "@inbook{Leon:1966,\n"
      "  title     = {Recent advances in optimization techniques: proceedings},\n"
      "  chapter   = {A comparison Among Eight Known Optimizing Procedures},\n"
      "  author    = {Leon, A.},\n"
      "  editor    = { Lavi, A. and Vogl, T.P.},\n"
      "  year      = {1966},\n"
      "  pages     = {28--46’,\n"
      "  publisher = {Wiley}\n"
      "}\n\n"
      "@article{doi:10.1137/0723046,\n"
      "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
      "  title   = {A Nonmonotone Line Search Technique for Newton’s Method},\n"
      "  journal = {SIAM Journal on Numerical Analysis},\n"
      "  year    = {1986},\n"
      "  volume  = {23},\n"
      "  number  = {4},\n"
      "  pages   = {707--716},\n"
      "  doi     = {10.1137/0723046},\n"
      "}\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & xx, int_type k ) const override {
    real_type x  = xx(0);
    real_type y  = xx(1);
    real_type x3 = x*x*x;
    switch ( k ) {
      case 0: return ((600*(x3-y)*x)+2)*x-2;
      case 1: return 200*(y-x3);
    }
    return 0;
  }

  void
  evalF( dvec_t const & xx, dvec_t & f ) const override {
    real_type x  = xx(0);
    real_type y  = xx(1);
    real_type x3 = x*x*x;
    f(0) = ((600*(x3-y)*x)+2)*x-2;
    f(1) = 200*(y-x3);
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
  jacobian( dvec_t const & xx, dvec_t & jac ) const override {
    real_type x  = xx(0);
    real_type y  = xx(1);
    real_type x2 = x*x;
    real_type x3 = x2*x;
    jac(0) = (3000*x3-1200*y)*x+2;
    jac(1) = jac(2) = -600*x2;
    jac(3) = 200;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = x(1) = 1;
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -1.2;
    x(1) = -1;
  }

  int_type
  numInitialPoint( ) const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
