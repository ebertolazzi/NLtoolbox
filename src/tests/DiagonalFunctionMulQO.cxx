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

class DiagonalFunctionMulQO : public nonlinearSystem {
public:

  DiagonalFunctionMulQO( int_type neq)
  : nonlinearSystem(
      "Diagonal Functions Multiplied by quasi-orthogonal matrix",
      "@article{Gasparo:2000,\n"
      "  Author    = {Maria Grazia Gasparo},\n"
      "  Title     = {A nonmonotone hybrid method for nonlinear systems},\n"
      "  Journal   = {Optimization Methods and Software},\n"
      "  Number    = {2},\n"
      "  Pages     = {79--94},\n"
      "  Publisher = {Taylor & Francis},\n"
      "  Volume    = {13},\n"
      "  Year      = {2000},\n"
      "  Doi       = {10.1080/10556780008805776},\n"
      "}\n",
      neq
    )
  { checkThree(n,3); }

  virtual
  real_type
  evalFk( dvec_t const & X, int_type i ) const override {
    int_type i3 = i/3;
    real_type x0 = X(i3*3+0);
    real_type x1 = X(i3*3+1);
    real_type x2 = X(i3*3+2);
    switch ( i % 3 ) {
      case 0: return x0*(0.6+1.6*x0*x0) + x1*(9.6-7.2*x1) - 4.8;
      case 1: return 0.48*x0+x1*(-4.32+x1*(3.24-0.72*x1))+x2*(0.2*x2*x2-1)+2.16;
      case 2: return x2*(1.25-0.25*x2*x2);
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & X, dvec_t & F ) const override {
    for ( int_type i = 0; i < n; i += 3 ) {
      real_type x0 = X(i+0);
      real_type x1 = X(i+1);
      real_type x2 = X(i+2);
      F(i+0) = x0*(0.6+1.6*x0*x0) + x1*(9.6-7.2*x1) - 4.8;
      F(i+1) = 0.48*x0+x1*(-4.32+x1*(3.24-0.72*x1))+x2*(0.2*x2*x2-1)+2.16;
      F(i+2) = x2*(1.25-0.25*x2*x2);
    }
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 2*n;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n; i += 3 ) {
      SETIJ(i+0,i+0);
      SETIJ(i+0,i+1);

      SETIJ(i+1,i+0);
      SETIJ(i+1,i+1);
      SETIJ(i+1,i+2);

      SETIJ(i+2,i+2);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; i += 3 ) {
      real_type x0 = X(i+0);
      real_type x1 = X(i+1);
      real_type x2 = X(i+2);

      jac(kk++) = 0.6+4.8*x0*x0;
      jac(kk++) = 9.6-14.4*x1;

      jac(kk++) = 0.48;
      jac(kk++) = -4.32+x1*(6.48-2.16*x1);
      jac(kk++) = 0.6*x2*x2-1;

      jac(kk++) = 1.25 -0.75*x2*x2;
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & X, int_type ) const override {
  }

  virtual
  void
  getInitialPoint( dvec_t & X, int_type ) const override {
    for ( int_type i = 0; i < n; i += 3 ) {
      X(i+0) = -1;
      X(i+1) = 0.5;
      X(i+2) = -1;
    }
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
