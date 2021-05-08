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

class CliffFunction : public nonlinearSystem {
public:

  CliffFunction()
  : nonlinearSystem(
      "Cliff Function",
      "@article{Grippo:1991,\n"
      "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
      "  title   = {A Class of Nonmonotone Stabilization Methods\n"
      "             in Unconstrained Optimization},\n"
      "  journal = {Numer. Math.},\n"
      "  year    = {1991},\n"
      "  volume  = {59},\n"
      "  number  = {1},\n"
      "  pages   = {779--805},\n"
      "  doi     = {10.1007/BF01385810},\n"
      "}\n",
      2
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type tmp = 1-20.0*exp(20.0*(x(0)-x(1)));
    switch ( k ) {
      case 0: return x(0)/5000.0-3.0/5000.0-tmp;
      case 1: return tmp;
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type tmp = 1-20.0*exp(20.0*(x(0)-x(1)));
    f(0) = x(0)/5000.0-3.0/5000.0-tmp;
    f(1) = tmp;
  }

  virtual
  int_type
  jacobianNnz() const override
  { return 4; }

  virtual
  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type tmp = 400*exp(20.0*(x(0)-x(1)));
    jac(0) = 1/5000.0+tmp;
    jac(1) = -tmp;
    jac(2) = -tmp;
    jac(3) = tmp;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 3;
    x(1) = 3 + log(20.0)/20.0;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = -1;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
