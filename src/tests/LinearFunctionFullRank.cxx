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

class LinearFunctionFullRank : public nonlinearSystem {

public:

  LinearFunctionFullRank()
  : nonlinearSystem(
      "Linear function - full rank",
      "@article{More:1981,\n"
      "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title   = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1981},\n"
      "  volume  = {7},\n"
      "  number  = {1},\n"
      "  pages   = {17--41},\n"
      "  doi     = {10.1145/355934.355936},\n"
      "}\n",
      10
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type sumx = 0;
    for ( int_type i = 0; i < n; ++i ) sumx += x(i);
    return x(k) - (2*sumx/n) - 1;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type sumx = 0;
    for ( int_type i = 0; i < n; ++i ) sumx += x(i);
    for ( int_type i = 0; i < n; ++i ) f(i) = x(i) - (2*sumx/n) - 1;
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
  jacobian( dvec_t const &, dvec_t & jac ) const override {
    real_type bf = real_type(2)/n;
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = 0; j < n; ++j ) {
        jac(kk) = -bf;
        if ( i == j ) jac(kk) += 1;
        ++kk;
      }
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.fill(1);
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = -1;
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
