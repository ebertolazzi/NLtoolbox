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

class Function27 : public nonlinearSystem {
public:

  Function27( integer neq)
  : nonlinearSystem(
      "Function 27",
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
    real_type f = 0;
    if ( i == 0 ) {
      f = x(0)*x(0);
      for (  i = 1; i < n; ++i ) f += x(i)*x(i);
    } else {
      f = -2*x(0)*x(i);
    }
    return f;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0)*x(0);
    for ( integer i = 1; i < n; ++i ) {
      f(0) += x(i)*x(i);
      f(i) = -2*x(0)*x(i);
    }
  }

  integer
  jacobianNnz() const override {
    return 3*n-2;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    for ( integer i = 1; i < n; ++i ) {
      SETIJ(0,i);
      SETIJ(i,0);
      SETIJ(i,i);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    jac(kk++) = 2*x(0);
    for ( integer i = 1; i < n; ++i ) {
      jac(kk++) = 2*x(i);
      jac(kk++) = -2*x(i);
      jac(kk++) = -2*x(0);
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
    x.fill(1.0/(n*n));
    x(0) = 100;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
