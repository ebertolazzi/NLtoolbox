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

class SingularFunction : public nonlinearSystem {

public:

  SingularFunction( int_type neq )
  : nonlinearSystem(
      "Singular Function",
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
  { checkMinEquations(n,2); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    if ( i == 0   ) return power3(x(0))/3+power2(x(1))/2;
    if ( i == n-1 ) return -power2(x(n-1))/2+n*power3(x(n-1))/3;
    return -x(i)*x(i)/2 + (i+1)*power3(x(i))/3 + power2(x(i+1))/2;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0)   = power3(x(0))/3+power2(x(1))/2;
    f(n-1) = power2(x(n-1))*((n/3.0)*x(n-1)-0.5);
    for ( int_type i = 1; i < n-1; ++i )
      f(i) = power2(x(i))*( ((i+1)/3.0)*x(i) - 0.5 ) + 0.5*power2(x(i+1));
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 2*(n-2)+3;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(n-1,n-1);
    for ( int_type i = 1; i < n-1; ++i ) {
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = power2(x(0));
    jac(kk++) = x(1);
    jac(kk++) = (n*x(n-1)-1)*x(n-1);
    for ( int_type i = 1; i < n-1; ++i ) {
      jac(kk++) = ((i+1)*x(i)-1)*x(i);
      jac(kk++) = x(i+1);
    }
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
    x.fill(1);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
      //ASSERT( std::abs(x(i)) < 10, "Bad range" );
      //ASSERT( x(i) >= 0 && x(i) < 1000, "Bad range" );
  }


};
