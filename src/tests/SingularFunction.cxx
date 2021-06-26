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

  SingularFunction( integer neq )
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

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( i == 0   ) return power3(x(0))/3+power2(x(1))/2;
    if ( i == n-1 ) return -power2(x(n-1))/2+n*power3(x(n-1))/3;
    return -x(i)*x(i)/2 + (i+1)*power3(x(i))/3 + power2(x(i+1))/2;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0)   = power3(x(0))/3+power2(x(1))/2;
    f(n-1) = power2(x(n-1))*((n/3.0)*x(n-1)-0.5);
    for ( integer i = 1; i < n-1; ++i )
      f(i) = power2(x(i))*( ((i+1)/3.0)*x(i) - 0.5 ) + 0.5*power2(x(i+1));
  }

  integer
  jacobianNnz() const override {
    return 2*(n-2)+3;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(n-1,n-1);
    for ( integer i = 1; i < n-1; ++i ) {
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    jac(kk++) = power2(x(0));
    jac(kk++) = x(1);
    jac(kk++) = (n*x(n-1)-1)*x(n-1);
    for ( integer i = 1; i < n-1; ++i ) {
      jac(kk++) = ((i+1)*x(i)-1)*x(i);
      jac(kk++) = x(i+1);
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.fill(1);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
      //ASSERT( abs(x(i)) < 10, "Bad range" );
      //ASSERT( x(i) >= 0 && x(i) < 1000, "Bad range" );
  }

};
