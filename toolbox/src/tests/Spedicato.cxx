/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

// AZZZ QUALCOSA NON TORNA @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SpedicatoFunction17 : public nonlinearSystem {
public:

  SpedicatoFunction17( int_type neq )
  : nonlinearSystem(
      "Spedicato N.17",
      "@Article{Spedicato1997,\n"
      "  author  = {Spedicato, E. and Huang, Z.},\n"
      "  title   = {Numerical experience with newton-like methods\n"
      "             for nonlinear algebraic systems},\n"
      "  journal = {Computing},\n"
      "  year    = {1997},\n"
      "  volume  = {58},\n"
      "  number  = {1},\n"
      "  pages   = {69--89},\n"
      "  doi     = {10.1007/BF02684472},\n"
      "}\n\n"
      "@book{meresoo:1990,\n"
      "  title     = {Test Examples of Systems of Nonlinear Equations: Version 3-90},\n"
      "  author    = {Meresoo, T. and Roose, A. and Kulla,\n"
      "               V. and Estonian Software and Computer Service Company},\n"
      "  year      = 1990,\n"
      "  publisher = {Estonian Software and Computer Service Company}\n"
      "}\n",
      neq
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    if ( i == 0   ) return x(0);
    if ( i == n-1 ) return x(n-1)-20;
    return x(i+1)+x(i)+x(i-1)+power2(x(i+1)-x(i-1))/4;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0)   = x(0);
    f(n-1) = x(n-1) - 20;
    for ( int_type i = 1; i < n-1; ++i )
      f(i) = x(i+1)+x(i)+x(i-1)+power2(x(i+1)-x(i-1))/4;
  }

  int_type
  jacobianNnz() const override {
    return 3*(n-2)+2;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(n-1,n-1);
    for ( int_type i = 1; i < n-1; ++i ) {
      SETIJ(i,i);
      SETIJ(i,i+1);
      SETIJ(i,i-1);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = 1;
    jac(kk++) = 1;
    for ( int_type i = 1; i < n-1; ++i ) {
      real_type bf = x(i+1)-x(i-1);
      jac(kk++) = 1;
      jac(kk++) = 1+bf/2;
      jac(kk++) = 1-bf/2;
    }
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 10;
  }

  int_type
  numInitialPoint() const override { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( int_type i = 0; i < n; ++i )
      NONLIN_ASSERT( abs(x(i)) < 10000, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(10000);
    L.fill(-10000);
  }

};
