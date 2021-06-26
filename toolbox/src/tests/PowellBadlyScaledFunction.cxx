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

class PowellBadlyScaledFunction : public nonlinearSystem {
  real_type const x0e;
  real_type const x1e;
  real_type const scale;
  
public:

  PowellBadlyScaledFunction()
  : nonlinearSystem(
      "Powell badly scaled function",
      "@inbook{Powell:1970,\n"
      "  title     = {Numerical methods for nonlinear algebraic equations},\n"
      "  booktitle = {Proceedings of a {C}onference, {U}niversity of {E}ssex,\n"
      "              {C}olchester, 6--7 {J}anuary 1969},\n"
      "  chapter   = {An hybrid method for non linear equations},\n"
      "  editor    = {Rabinowitz, Philip},\n"
      "  publisher = {Gordon and Breach Science Publishers, London-New York-Paris},\n"
      "  year      = {1970},\n"
      "  pages     = {87--114},\n"
      "}\n\n"
      "@article{More:1981,\n"
      "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title   = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1981},\n"
      "  volume  = {7},\n"
      "  number  = {1},\n"
      "  pages   = {17--41},\n"
      "  doi     = {10.1145/355934.355936},\n"
      "}\n",
      2
    )
  //, x0e(1e-6)
  //, x1e(1e6)
  //, scale(1)
  //, x0e(0.109815932969981745568376164563E-4)
  //, x1e(9.10614673986652401094671049032e8)
  //, scale(10000e-18)
  , x0e(0.109815932969981745568376164563E-4)
  , x1e(9.10614673986652401094671049032)
  , scale(10000)
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = scale*((x(0) * x(1)) - (x0e*x1e));
    f(1) = (exp(-x(1))-exp(-x1e)) + (exp(-x(0))-exp(-x0e));
    //f(0) = 10000*(x(0) * x(1)) - 1;
    //f(1) = exp(-x(1)) + exp(-x(0)) - 1.0001;
  }

  integer
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
    jac(0) = scale * x(1);
    jac(1) = scale * x(0);

    jac(2) = -exp( -x(0) );
    jac(3) = -exp( -x(1) );
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    //x(0) = 0.109815932969981745568376164563E-4;
    //x(1) = 9.10614673986652401094671049032;
    x(0) = x0e;
    x(1) = x1e;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 100;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  UTILS_ASSERT(abs(x(i)) < 1000, "Bad range" );
    // UTILS_ASSERT( x(i) >= 0i] < 1000, "Bad range" );
  }

};
