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

class ZeroJacobianFunction : public nonlinearSystem {
public:

  ZeroJacobianFunction( int_type neq )
  : nonlinearSystem(
      "Zero Jacobian Function (same as function 27)",
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
  {  checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type f = 0;
    if ( i == 0 ) {
      f += x.dot(x);
    } else {
      f = -2*x(0)*x(i);
    }
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x.dot(x);
    for ( int_type i = 1; i < n; ++i ) f(i)  = -2*x(0)*x(i);
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
    for ( int_type i = 0; i < n; ++i ) { SETIJ(0,i); }
    for ( int_type i = 1; i < n; ++i ) {
      SETIJ(i,0);
      SETIJ(i,i);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) jac(kk++) = 2*x(i);
    for ( int_type i = 1; i < n; ++i ) {
      jac(kk++) = -2*x(i);
      jac(kk++) = -2*x(0);
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override
  { }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    real_type bf = (1.0/60.0-100.0/(6.0*n))*(1.0/60.0-50.0/6.0);
    x.fill(bf);
    x(0) = 100.0*(n-100.0)/n;
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