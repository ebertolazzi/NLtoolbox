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

class CompressibilityFactorFromTheRKequation : public nonlinearSystem {
  real_type Q, r;
public:

  CompressibilityFactorFromTheRKequation()
  : nonlinearSystem(
      "Compressibility factor from the RK equation",
      "@book{Cutlip:2007,\n"
      "  author    = {Cutlip, Michael and Shacham, Mordechai},\n"
      "  title     = {Problem Solving in Chemical and Biochemical Engineering\n"
      "               with Polymath,\\texttrademark Excel,\n"
      "               and Matlab\\textregistered, Second Edition},\n"
      "  year      = {2007},\n"
      "  isbn      = {9780131482043},\n"
      "  publisher = {Prentice Hall Press},\n"
      "}\n",
      1
    )
  {
    real_type P    = 200;
    real_type Pc   = 33.5;
    real_type T    = 631*2;
    real_type Tc   = 126.2;
    real_type Pr   = P / Pc;
    real_type Tr   = T / Tc;
    real_type Asqr = 0.4278*Pr/pow(Tr,2.5);
    real_type B    = 0.0867*Pr/Tr;
	  Q = B*B+B-Asqr;
    r = Asqr*B;
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    real_type z = x(0);
    return ((z-1)*z - Q)*z - r;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type z = x(0);
    f(0) = ((z-1)*z - Q)*z - r;
  }

  integer
  jacobianNnz() const override
  { return 1; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type z = x(0);
    jac(0) = (3*z-2)*z - Q;
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0.65;  break;
      case 1: x(0) = -0.5;  break;
      case 2: x(0) = 1;     break;
      case 3: x(0) = -0.02; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 4; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
