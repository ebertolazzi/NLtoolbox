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

class BrownAlmostLinearFunction : public nonlinearSystem {
public:
  
  BrownAlmostLinearFunction( int_type  neq )
  : nonlinearSystem(
      "Brown almost linear function",
      "@article{Brown:1968,\n"
      "  author    = {Brown, Kenneth M.},\n"
      "  title     = {A Quadratically Convergent Newton-Like Method Based\n"
      "               Upon Gaussian-Elimination},\n"
      "  volume    = {6},\n"
      "  number    = {4},\n"
      "  year      = {1969},\n"
      "  pages     = {560--569}\n"
      "  publisher = {SIAM Journal on Numerical Analysis},\n"
      "}\n\n"
      "@article{More:1981,\n"
      "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title   = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  volume  = {7},\n"
      "  number  = {1},\n"
      "  year    = {1981},\n"
      "  pages   = {17--41},\n"
      "  doi     = {10.1145/355934.355936},\n"
      "}\n",
      neq
    )
  { checkMinEquations(neq,2); }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    if ( k < n-1 ) {
      real_type sumx = 0;
      for ( int_type i = 0; i < n; ++i ) sumx  += x(i);
      return x(k) + (sumx - (n+1));
    } else {
      real_type prodx = 1;
      for ( int_type i = 0; i < n; ++i ) prodx *= x(i);
      return prodx - 1;
    }           
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type sumx  = x.sum();
    real_type prodx = x.prod();
    for ( int_type i = 0; i < n-1; ++i ) f(i) = x(i) + (sumx - (n+1));
    f(n-1) = prodx - 1;
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0; // fortran address
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {

    jac.fill(1);
    for ( int_type i = 0; i < n-1; ++i ) jac(caddr(i,i)) = 2;

    // last row
    for ( int_type j = 0; j < n; ++j ) {
      real_type prod = 1;
      for ( int_type k = 0; k < n; ++k ) {
        if ( k != j ) prod *= x(k);
      }
      jac(caddr(n-1,j)) = prod;
    }
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.fill(1);
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill( 0.5 );
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};

