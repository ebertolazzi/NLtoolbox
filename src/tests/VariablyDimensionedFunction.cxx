/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

// TEST 221
/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class VariablyDimensionedFunction : public nonlinearSystem {
  mutable dvec_t gf2;
public:
  
  VariablyDimensionedFunction( int_type neq )
  : nonlinearSystem(
      "Variably dimensioned function",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n"
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
      neq
    )
  { checkMinEquations(neq,2); gf2.resize(n); }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type sum1 = 0;
    for ( int_type j = 0; j < n; ++j )
      sum1 += (j+1)*(x(j)-1);
    return x(k) - 1 + (k+1)*sum1*(1+2*power2(sum1));
  }

  //  f = f1 * f1 * ( 1.0 + f1 * f1 ) + f2;

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type sum1 = 0;
    for ( int_type j = 0; j < n; ++j )
      sum1 += (j+1)*(x(j)-1);
    for ( int_type j = 0; j < n; ++j )
      f(j) = x(j) - 1 + (j+1)*sum1*(1+2*power2(sum1));
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type sum1 = 0;
    for ( int_type j = 0; j < n; ++j )
      sum1 += (j+1)*(x(j)-1);

    int_type kk = 0;
    for ( int_type k = 0; k < n; ++k ) {
    	for ( int_type j = 0; j < n; ++j ) {
    	  jac(kk) = (k+1)*(1+6*power2(sum1))*(j+1);
        if ( j == k ) jac(kk) += 1;
        ++kk;
      }
    }
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }
  
  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k )
      x(k) = 1 - real_type(k+1) /real_type(n);
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
