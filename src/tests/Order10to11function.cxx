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

class Order10to11function : public nonlinearSystem {
  real_type const SCALE;
public:

  Order10to11function()
  : nonlinearSystem(
      "Order10to11 function",
      "@article{Shacham:1972,\n"
      "  author  = {Mordechai Shacham and Ephraim Kehat},\n"
      "  title   = {An iteration method with memory for\n"
      "             the solution of a non-linear equation},\n"
      "  journal = {Chemical Engineering Science},\n"
      "  volume  = {27},\n"
      "  number  = {11},\n"
      "  pages   = {2099--2101},\n"
      "  year    = {1972},\n"
      "  doi     = {10.1016/0009-2509(72)87067-2}\n"
      "}\n",
      1
    )
  , SCALE(1e-10)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type T = x(0);
    return SCALE*(exp(21000./T)/(T*T) - 1.11E11);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type T = x(0);
    f(0) = SCALE*(exp(21000./T)/(T*T) - 1.11E11);
  }

  int_type
  jacobianNnz() const override
  { return 1; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type T = x(0);
    real_type tmp = 21000/T;
    jac(0) = -SCALE*(exp(tmp)/power3(T))*(2+tmp);
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 555;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
