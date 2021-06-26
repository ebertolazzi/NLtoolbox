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

#define EXPONENTIAL_FUNCTION_BIBTEX \
"@article{LaCruz:2003,\n" \
"  author    = { William {La Cruz}  and  Marcos Raydan},\n" \
"  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems},\n" \
"  journal   = {Optimization Methods and Software},\n" \
"  year      = {2003},\n" \
"  volume    = {18},\n" \
"  number    = {5},\n" \
"  pages     = {583--599},\n" \
"  publisher = {Taylor & Francis},\n" \
"  doi       = {10.1080/10556780310001610493},\n" \
"}\n"

class ExponentialFunction1 : public nonlinearSystem {
public:

  ExponentialFunction1( integer neq )
  : nonlinearSystem( "Exponential Function N.1", EXPONENTIAL_FUNCTION_BIBTEX, neq )
  { checkMinEquations(n,1); }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( i == 0 ) return exp(x(0)-1) - 1;
    return (i+1)*(exp(x(i)-1)-x(i));
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = exp(x(0)-1) - 1;
    for ( integer  i = 1; i < n; ++i )
      f(i) = (i+1)*(exp(x(i)-1)-x(i));
  }

  integer
  jacobianNnz() const override
  { return n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      { ii(kk) = jj(kk) = i; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = exp(x(0)-1);
    for ( integer i = 1; i < n; ++i )
      jac(i) = (i+1)*(exp(x(i)-1)-1);
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override
  { }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.fill( n/(n-1.0) );
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ExponentialFunction2 : public nonlinearSystem {
public:

  ExponentialFunction2( integer neq )
  : nonlinearSystem( "Exponential Function N.2", EXPONENTIAL_FUNCTION_BIBTEX, neq )
  { checkMinEquations(n,1); }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( i == 0 ) return exp(x(0)) - 1;
    return ((i+1)/10.0)*(exp(x(i))+x(i-1)-1);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const  override{
    f(0) = exp(x(0)) - 1;
    for ( integer i = 1; i < n; ++i )
      f(i) = ((i+1)/10.0)*(exp(x(i))+x(i-1)-1);
  }

  integer
  jacobianNnz() const override
  { return 2*n-1; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    ii(kk) = jj(kk) = 0; ++kk;
    for ( integer i = 1; i < n; ++i ) {
      ii(kk) = jj(kk) = i; ++kk;
      ii(kk) = i; jj(kk) = i-1; ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    jac(kk++) = exp(x(0));
    for ( integer i = 1; i < n; ++i ) {
      jac(kk++) = ((i+1)/10.0)*exp(x(i));
      jac(kk++) = ((i+1)/10.0);
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
    x.fill( 1.0/(n*n) );
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class ExponentialFunction3 : public nonlinearSystem {
public:

  ExponentialFunction3( integer neq )
  : nonlinearSystem( "Exponential Function N.3", EXPONENTIAL_FUNCTION_BIBTEX, neq )
  { checkMinEquations(n,1); }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( i == n-1 ) return (0.1*n)*(1-exp(-x(n-1)*x(n-1)));
    return (0.1*(i+1))*(1-x(i)*x(i)-exp(-x(i)*x(i)));
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(n-1) = (0.1*n)*(1-exp(-x(n-1)*x(n-1)));
    for ( integer i = 0; i < n-1; ++i )
      f(i) = (0.1*(i+1))*(1-x(i)*x(i)-exp(-x(i)*x(i)));
  }

  integer
  jacobianNnz() const override
  { return n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      { ii(kk) = jj(kk) = i; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n-1; ++i )
      jac(kk++) = 0.2*(i+1)*x(i)*(exp(-x(i)*x(i))-1);
    jac(kk++) = 0.2*n*x(n-1)*exp(-x(n-1)*x(n-1));
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    real_type bf = 1.0/(4.0*n*n);
    for ( integer i = 0; i < n; ++i ) x(i) = (i+1)*bf;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
