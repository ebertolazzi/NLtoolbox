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

class LogarithmicFunction : public nonlinearSystem {
public:

  LogarithmicFunction( int_type neq)
  : nonlinearSystem(
      "Logarithmic Function",
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
  { checkMinEquations(n,1); }

  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    return log(x(i)+1)-x(i)/n;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n; ++i )
      f(i) = log(x(i)+1)-x(i)/n;
  }

  int_type
  jacobianNnz() const override
  { return n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    for ( int_type i = 0; i < n; ++i ) ii(i) = jj(i) = i;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    for ( int_type i = 0; i < n; ++i )
      jac(i) = 1.0/(x(i)+1)-1.0/n;
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override
  { }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = n+1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
