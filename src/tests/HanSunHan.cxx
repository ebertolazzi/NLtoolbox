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

class HanSunHan : public nonlinearSystem {
  real_type rr[99];
public:

  HanSunHan()
  : nonlinearSystem(
      "Han-Sun-Han-SAMPAJO 2005 function test",
      "@article{Han:2005,\n"
      "  author  = {Qiaoming Han and Wenyu Sun and Jiye Han and Raimudo J. B. Sampaio},\n"
      "  title   = {An adaptive conic trust-region method for unconstrained optimization},\n"
      "  journal = {Optimization Methods and Software},\n"
      "  year    = {2005},\n"
      "  volume  = {20},\n"
      "  number  = {6},\n"
      "  pages   = {665--677},\n"
      "  doi = {10.1080/10556780410001697677}\n"
      "}\n",
      2
    )
  {
    for ( integer i = 0; i < 99; ++i ) {
      real_type arg = (i+1) / 100.0;
      rr[i] = pow(-50.0*log(arg),2.0/3.0) + 25.0;
    }
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type t1  = x(0)*x(0);
    real_type t3  = t1*t1;
    real_type t7  = x(1)*x(1);
    real_type t12 = exp(t7);
    f(0) = (8 * t3 * t1 * x(0)) + (2 * x(0)) + 2 * x(0) * t7;
    f(1) = (2 * t1 * x(1)) + 2 * x(1) * t12;
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
    real_type t1  = x(0)*x(0);
    real_type t2  = t1*t1;
    real_type t5  = x(1)*x(1);
    real_type t9  = 4*x(0)*x(1);
    real_type t11 = exp(t5);
    jac(0) = (56 * t2 * t1) + 2 + 2 * t5;
    jac(1) = jac(2) = t9;
    jac(3) = (2 * t1) + 2 * t11 + 4 * t5 * t11;
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 5;
    x(1) = 3;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
