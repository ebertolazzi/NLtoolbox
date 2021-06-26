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

class Function21 : public nonlinearSystem {
public:

  Function21( integer neq)
  : nonlinearSystem(
      "Function 21",
      "@article{LaCruz:2003,\n"
      "  author    = {William {La Cruz}  and  Marcos Raydan},\n"
      "  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems},\n"
      "  journal   = {Optimization Methods and Software},\n"
      "  year      = {2003},\n"
      "  volume    = {18},\n"
      "  number    = {5},\n"
      "  pages     = {583--599},\n"
      "  publisher = {Taylor & Francis},\n"
      "  doi       = {10.1080/10556780310001610493},\n"
      "}\n",
      neq
    )
  { checkThree(n,3); }

  real_type
  evalFk( dvec_t const & X, integer i ) const override {
    integer i3 = i/3;
    real_type const & x = X(i3*3+0);
    real_type const & y = X(i3*3+1);
    real_type const & z = X(i3*3+2);
    switch ( i % 3 ) {
      case 0: return x*y - z*z - 1;
      case 1: return x*y*z - x*x + y*y - 2;
      case 2: return exp(-x)-exp(-y);
    }
    return 0;
  }

  void
  evalF( dvec_t const & X, dvec_t & F ) const override {
    for ( integer i = 0; i < n; i += 3 ) {
      real_type const & x = X(i+0);
      real_type const & y = X(i+1);
      real_type const & z = X(i+2);
      F(i+0) = x*y - z*z - 1;
      F(i+1) = x*y*z - x*x + y*y - 2;
      F(i+2) = exp(-x)-exp(-y);
    }
  }

  integer
  jacobianNnz() const override {
    return 3*n;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( integer i = 0; i < n; i += 3 ) {
      integer I0 = i+0;
      integer I1 = i+1;
      integer I2 = i+2;

      SETIJ(I0,I0);
      SETIJ(I0,I1);
      SETIJ(I0,I2);

      SETIJ(I1,I0);
      SETIJ(I1,I1);
      SETIJ(I1,I2);

      SETIJ(I2,I0);
      SETIJ(I2,I1);
      SETIJ(I2,I2);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; i += 3 ) {

      integer I0 = i+0;
      integer I1 = i+1;
      integer I2 = i+2;

      real_type const & x = X(I0);
      real_type const & y = X(I1);
      real_type const & z = X(I2);

      jac(kk++) = y;
      jac(kk++) = x;
      jac(kk++) = -2*z;

      jac(kk++) = y*z - 2*x;
      jac(kk++) = x*z + 2*y;
      jac(kk++) = x*y;

      jac(kk++) = -exp(-x);
      jac(kk++) = exp(-y);
      jac(kk++) = 0;
    }
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    for ( integer i = 0; i < n; ++i ) x(i) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
