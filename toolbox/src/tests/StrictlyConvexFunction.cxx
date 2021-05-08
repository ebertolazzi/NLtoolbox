/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define STRICT_CONVEX_FUNCTION_BIBTEX \
"@article{Raydan:1997,\n" \
"  author  = {Raydan, M.},\n" \
"  title   = {The Barzilai and Borwein Gradient Method for\n" \
"             the Large Scale Unconstrained Minimization Problem},\n" \
"  journal = {SIAM Journal on Optimization},\n" \
"  volume  = {7},\n" \
"  number  = {1},\n" \
"  pages   = {26-33},\n" \
"  year    = {1997},\n" \
"  doi     = {10.1137/S1052623494266365},\n" \
"}\n\n" \
"@article{LaCruz:2003,\n" \
"  author    = {William {La Cruz}  and  Marcos Raydan},\n" \
"  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems},\n" \
"  journal   = {Optimization Methods and Software},\n" \
"  year      = {2003},\n" \
"  volume    = {18},\n" \
"  number    = {5},\n" \
"  pages     = {583--599},\n" \
"  publisher = {Taylor & Francis},\n" \
"  doi       = {10.1080/10556780310001610493},\n" \
"}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class StrictlyConvexFunction1 : public nonlinearSystem {
public:

  StrictlyConvexFunction1( int_type neq )
  : nonlinearSystem(
      "Strictly Convex Function 1",
      STRICT_CONVEX_FUNCTION_BIBTEX,
      neq
    )
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    return exp(x(i))-1;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n; ++i )
      f(i) = exp(x(i))-1;
  }

  virtual
  int_type
  jacobianNnz() const override {
    return n;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n; ++i ) { SETIJ(i,i); }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    for ( int_type i = 0; i < n; ++i ) jac(i) = exp(x(i));
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = (i+1.0)/n;
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class StrictlyConvexFunction2 : public nonlinearSystem {
public:

  StrictlyConvexFunction2( int_type neq )
  : nonlinearSystem(
      "Strictly Convex Function 2",
      STRICT_CONVEX_FUNCTION_BIBTEX,
      neq
    )
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    return ((i+1.0)/10.0)*(exp(x(i))-1);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n; ++i )
      f(i) = ((i+1.0)/10.0)*(exp(x(i))-1);
  }

  virtual
  int_type
  jacobianNnz() const override {
    return n;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n; ++i ) { SETIJ(i,i); }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    for ( int_type i = 0; i < n; ++i )
      jac(i) = ((i+1.0)/10.0)*exp(x(i));
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 1;
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
