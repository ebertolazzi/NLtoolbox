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

class Function18 : public nonlinearSystem {
public:

  Function18( integer neq)
  : nonlinearSystem(
      "Function 18",
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
    dvec_t const & x = X.segment(i3*3,3);
    switch ( i % 3 ) {
      case 0: return x(0)*x(1) - x(2)*x(2) - 1;
      case 1: return x(0)*x(1)*x(2) - x(0)*x(0) + x(1)*x(1) - 2;
      case 2: return exp(x(0))-exp(x(1));
    }
    return 0;
  }

  void
  evalF( dvec_t const & X, dvec_t & F ) const override {
    for ( integer i = 0; i < n; i += 3 ) {
      dvec_t const & x = X.segment(i,3);
      F(i+0) = x(0)*x(1) - x(2)*x(2) - 1;
      F(i+1) = x(0)*x(1)*x(2) - x(0)*x(0) + x(1)*x(1) - 2;
      F(i+2) = exp(x(0))-exp(x(1));
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
      SETIJ(i+0,i+0);
      SETIJ(i+0,i+1);
      SETIJ(i+0,i+2);

      SETIJ(i+1,i+0);
      SETIJ(i+1,i+1);
      SETIJ(i+1,i+2);

      SETIJ(i+2,i+0);
      SETIJ(i+2,i+1);
      SETIJ(i+2,i+2);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; i += 3 ) {
      dvec_t const & x = X.segment(i,3);

      jac(kk++) = x(1);
      jac(kk++) = x(0);
      jac(kk++) = -2*x(2);

      jac(kk++) = x(1)*x(2) - 2*x(0);
      jac(kk++) = x(0)*x(2) + 2*x(1);
      jac(kk++) = x(0)*x(1);

      jac(kk++) = exp(x(0));
      jac(kk++) = -exp(x(1));
      jac(kk++) = 0;
    }
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & X, integer ) const override {
    for ( integer i = 0; i < n; i += 3 ) {
      X(i+0) = sqrt(2.0);
      X(i+1) = sqrt(2.0);
      X(i+2) = 1;
    }
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.setZero();
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
