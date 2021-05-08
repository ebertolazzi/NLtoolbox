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

class Hilbert : public nonlinearSystem {
public:

  Hilbert( int_type n )
  : nonlinearSystem(
      "Hilbert Matrix Function F = x'Ax",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      n
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type f = 0;
    for ( int_type j = 0; j < n; ++j )
      f += (2.0 * x(k)) / ( k + j + 1 );
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n; ++i ) {
      f(i) = 0;
      for ( int_type j = 0; j < n; ++j ) {
        f(i) += (2.0 * x(j)) / ( i + j + 1 );
      }
    }
  }

  virtual
  int_type
  jacobianNnz() const override
  { return n*n; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = 0; j < n; ++j ) {
        jac(kk++) = 2.0 / ( i + j + 1 );
      }
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.setZero();
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

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
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
