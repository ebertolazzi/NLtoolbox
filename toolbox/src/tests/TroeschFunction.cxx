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

class TroeschFunction : public nonlinearSystem {

  real_type const rho;
  real_type const h;

public:

  TroeschFunction( int_type neq )
  : nonlinearSystem(
      "Troesch Function",
      "@inproceedings{Varadhan2009,\n"
      "  author={R. Varadhan and Paul D. Gilbert},\n"
      "  title={{BB:} An {R} Package for Solving a Large System of\n"
      "         Nonlinear Equations and for Optimizing a High-Dimensional\n"
      "         Nonlinear Objective Function},\n"
      "  year={2009}\n"
      "}\n",
      neq
    )
  , rho(10)
  , h(1.0/(neq+1))
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type bf = rho*h*h;
    real_type f  = 2*x(i) + bf*sinh(rho*x(i));
    if      ( i == 0   ) f -= x(1);
    else if ( i == n-1 ) f -= x(n-2)+1;
    else                 f -= x(i-1) + x(i+1);
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type bf = rho*h*h;

    for ( int_type i = 0; i < n; ++i )
      f(i) = 2*x(i) + bf*sinh(rho*x(i));

    f(0)   -= x(1);
    f(n-1) -= x(n-2)+1;

    for ( int_type i = 1; i < n-1; ++i )
      f(i) -= x(i-1) + x(i+1);
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n;   ++i ) { SETIJ(i,i); }
    for ( int_type i = 0; i < n-1; ++i ) { SETIJ(i,i+1); }
    for ( int_type i = 1; i < n;   ++i ) { SETIJ(i,i-1); }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    real_type bf = rho*rho*h*h;
    for ( int_type i = 0; i < n;   ++i ) jac(kk++) = 2 + bf*cosh(rho*x(i));
    for ( int_type i = 0; i < n-1; ++i ) jac(kk++) = -1;
    for ( int_type i = 1; i < n;   ++i ) jac(kk++) = -1;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i )
      x(i) = ((123*i)%1001)/1000.0;
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
