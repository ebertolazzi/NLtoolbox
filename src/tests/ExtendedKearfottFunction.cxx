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

class ExtendedKearfottFunction : public nonlinearSystem {
public:

  ExtendedKearfottFunction()
  : nonlinearSystem(
      "Extended Kearfott Function",
      "@phdthesis{Kearfott:1977,\n"
      "  author    = {Kearfott, Ralph Baker},\n"
      "  title     = {Computing the Degree of Maps and a Generalized\n"
      "               Method of Bisection},\n"
      "  year      = {1977},\n"
      "  note      = {AAI7723103},\n"
      "  publisher = {The University of Utah},\n"
      "}\n\n"
      "@Article{Kearfott:1979,\n"
      "  author  = {Kearfott, Baker},\n"
      "  title   = {An efficient degree-computation method for a generalized\n"
      "             method of bisection},\n"
      "  journal = {Numerische Mathematik},\n"
      "  year    = {1979},\n"
      "  volume  = {32},\n"
      "  number  = {2},\n"
      "  pages   = {109--127},\n"
      "  doi     = {10.1007/BF01404868}\n"
      "}\n",
      7
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    if ( i == n-1 ) return x(n-1)*x(n-1) - x(0);
    return x(i)*x(i) - x(i+1);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n-1; ++i )
      f(i) = x(i)*x(i) - x(i+1);
    f(n-1) = x(n-1)*x(n-1) - x(0);
  }

  virtual
  int_type
  jacobianNnz() const override
  { return 2*n; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n-1; ++i ) {
      ii(kk) = jj(kk) = i; ++kk;
      ii(kk) = i; jj(kk) = i+1; ++kk;
    }
    ii(kk) = jj(kk) = n-1; ++kk;
    ii(kk) = n-1; jj(kk) = 0;
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n-1; ++i ) {
      jac(kk++) = 2*x(i);
      jac(kk++) = -1;
    }
    jac(kk++) = 2*x(n-1);
    jac(kk++) = -1;
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
    x.fill(0.1);
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
