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

class Colville : public nonlinearSystem {
public:

  Colville()
  : nonlinearSystem(
      "Colville Polynomial",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      4
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);

    real_type x1_2 = x1*x1;
    real_type x1_3 = x1_2*x1;
    real_type x3_2 = x3*x3;
    real_type x3_3 = x3_2*x3;
    
    f(0) =  400.0 * x1_3  - 400.0 * x2 * x1 + 2.0  * x1 - 2.0;
    f(1) = -200.0 * x1_2  + 220.2 * x2      + 19.8 * x4 - 40.0;
    f(2) = -360.0 * x3*x4 + 360.0 * x3_3    + 2.0  * x3 - 2.0;
    f(3) =  180.0 * x4    - 180.0 * x3_2    + 20.2 * x4 + 19.8 * x2 - 40.0;

  }

  virtual
  int_type
  jacobianNnz() const override
  { return 10; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I-1; jj(kk) = J-1; ++kk

    SETIJ(1,1);
    SETIJ(1,2);

    SETIJ(2,1);
    SETIJ(2,2);
    SETIJ(2,4);

    SETIJ(3,3);
    SETIJ(3,4);

    SETIJ(4,2);
    SETIJ(4,3);
    SETIJ(4,4);

    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);

    jac(0) = 1200.0 * x1*x1 - 400.0 * x2 + 2.0;
    jac(1) = -400.0 * x1;

    jac(2) = -400.0 * x1;
    jac(3) = 220.2;
    jac(4) = 19.8;

    jac(5) = -360.0 * x4 + 1080.0 * x3*x3 + 2.0;
    jac(6) = -360.0 * x3;

    jac(7) = 19.8;
    jac(8) = -360.0 * x3;
    jac(9) = 200.2;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 1;
    x(2) = 1;
    x(3) = 1;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.5;
    x(1) = 1.0;
    x(2) = -0.5;
    x(3) = -1.0;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( std::abs(x(i)) < 10000, "x[" << i << "] = "<< x(i) << " too big");
  }

};
