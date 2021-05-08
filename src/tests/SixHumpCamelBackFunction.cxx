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

class SixHumpCamelBackFunction : public nonlinearSystem {

public:

  SixHumpCamelBackFunction()
  : nonlinearSystem(
      "Six Hump Camel Back function",
      "Molga M. and Smutnicki C. (2005).\n"
      "Test functions for optimization needs,\n"
      "http://www.zsd.ict.pwr.wroc.pl/files/docs/functions\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x_in, int_type k ) const override {
    real_type x   = x_in[0];
    real_type y   = x_in[1];
    real_type t2  = x*x;
    real_type t5  = t2*t2;
    real_type t10 = y*y;
    switch ( k ) {
      case 0: return 8.0*x-0.84E1*t2*x+2.0*t5*x+y;
      case 1: return x-8.0*y+16.0*t10*y;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x   = x_in[0];
    real_type y   = x_in[1];
    real_type t2  = x*x;
    real_type t5  = t2*t2;
    real_type t10 = y*y;
    f(0) = 8.0*x-0.84E1*t2*x+2.0*t5*x+y;
    f(1) = x-8.0*y+16.0*t10*y;
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
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x  = x_in[0];
    real_type y  = x_in[1];
    real_type t1 = x*x;
    real_type t3 = t1*t1;
    real_type t6 = y*y;
    jac(0) = 8.0-0.252E2*t1+10.0*t3;
    jac(1) = jac(2) = 1.0;
    jac(3) = -8.0+48.0*t6;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = -0.08984;
    x(1) = 0.71266;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -5;
    x(1) = 5;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
