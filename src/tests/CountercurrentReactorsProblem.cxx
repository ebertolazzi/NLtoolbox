/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define COUNTERCURRENT_BIBTEX \
"@article{Bogle:1990,\n" \
"  author  = {Bogle, I. and Perkins, J.},\n" \
"  title   = {A New Sparsity Preserving Quasi-Newton Update\n" \
"             for Solving Nonlinear Equations},\n" \
"  journal = {SIAM Journal on Scientific and Statistical Computing},\n" \
"  year    = {1990},\n" \
"  volume  = {11},\n" \
"  number  = {4},\n" \
"  pages   = {621-630},\n" \
"  doi     = {10.1137/0911036},\n" \
"}\n" \
"@techreport{Bodon:1990,\n" \
"  author  = {Elena Bodon and Ladislav Luksan and Emilio Spedicato},\n" \
"  title   = {Numerical performance of ABS codes for nonlinear least squares},\n" \
"  year    = {2001},\n" \
"  number  = {Tech. Rep. DMSIA 27/2001, Universita degli Studi di Bergamo}\n" \
"}\n"

//  Problem N.2

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class CountercurrentReactorsProblem1 : public nonlinearSystem {
  real_type const alpha;
  real_type const theta;
public:
 
  CountercurrentReactorsProblem1( int_type neq )
  : nonlinearSystem(
      "Countercurrent Reactors Problem N.1",
      COUNTERCURRENT_BIBTEX,
      neq
    )
  , alpha(0.5)
  , theta(4.0)
  { checkMinEquations(neq,4); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const {
    for ( int_type i = 0; i < n; i += 2 ) {
      real_type xm2, xm1, xp2, xp3;
      if ( i == 0 ) {
        xm2 = 1;
        xm1 = 0;
      } else {
        xm2 = x(i-2);
        xm1 = x(i-1);
      }
      if ( i >= n-2 ) {
         xp2 = 0;
         xp3 = 1;
      } else {
         xp2 = x(i+2);
         xp3 = x(i+3);
      }
      real_type xi  = x(i);
      real_type xp1 = x(i+1);
      f(i)   = alpha * xm2 + (alpha-1)*xp2 - xi*(1+theta*xp1);
      f(i+1) = (alpha-1) * xm1 + (alpha-2)*xp3 - theta*xi*xp1;
    }
  }

  virtual
  int_type
  jacobianNnz() const {
    int_type kk = 0;
    for ( int_type i = 0; i < n; i += 2 ) {
      if ( i > 0   ) kk += 2;
      if ( i < n-3 ) kk += 2;
      kk += 4;
    }
    return kk;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const {
    int_type kk = 0;
    for ( int_type i = 0; i < n; i += 2 ) {
      if ( i > 0 ) {
        ii(kk) = i;   jj(kk) = i-2; ++kk;
        ii(kk) = i+1; jj(kk) = i-1; ++kk;
      }
      if ( i < n-3 ) {
        ii(kk) = i;   jj(kk) = i+2; ++kk;
        ii(kk) = i+1; jj(kk) = i+3; ++kk;
      }
      ii(kk) = i;   jj(kk) = i;   ++kk;
      ii(kk) = i;   jj(kk) = i+1; ++kk;
      ii(kk) = i+1; jj(kk) = i;   ++kk;
      ii(kk) = i+1; jj(kk) = i+1; ++kk;
    }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const {
    int_type kk = 0;
    for ( int_type i = 0; i < n; i += 2 ) {
      if ( i > 0 ) {
        jac(kk++) = alpha;
        jac(kk++) = alpha-1;
      }
      if ( i < n-3 ) {
        jac(kk++) = alpha-1;
        jac(kk++) = alpha-2;
      }
      real_type xi  = x(i);
      real_type xp1 = x(i+1);
      jac(kk++) = -(1+theta*xp1);
      jac(kk++) = -theta*xi;
      jac(kk++) = -theta*xp1;
      jac(kk++) = -theta*xi;
    }
  }

  virtual
  int_type
  numExactSolution() const { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const {
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type istart ) const {
    for ( int_type i = 0; i < n; ++i ) {
      switch ( i % 8 ) {
        case 0:         x(i) = 0.1; break;
        case 1: case 7: x(i) = 0.2; break;
        case 2: case 6: x(i) = 0.3; break;
        case 3: case 5: x(i) = 0.4; break;
        case 4:         x(i) = 0.5; break;
      }
    }
    switch ( istart ) {
    case 1: for ( int_type i = 0; i < n; ++i ) x(i) *= 10;  break;
    case 2: for ( int_type i = 0; i < n; ++i ) x(i) *= 100; break;
    }
  }

  virtual
  int_type
  numInitialPoint() const
  { return 3; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( std::abs(x(i)) < 1000000, "x range" );
  }

};


class CountercurrentReactorsProblem2 : public nonlinearSystem {
  real_type const A0;
  real_type const A1;
  real_type const B0;
  real_type const theta;
public:
 
  CountercurrentReactorsProblem2( int_type neq )
  : nonlinearSystem(
      "Countercurrent Reactors Problem N.2",
      COUNTERCURRENT_BIBTEX,
      neq
    )
  , A0(1)
  , A1(0.414214)
  , B0(0)
  , theta(4)
  { checkMinEquations(neq,6); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = A0*x(0) - (1-x(0))*x(2) - A1 - theta*A1*x(1);
    f(1) = B0*x(0) - (1-x(0))*x(3) - A1 - theta*A1*x(1);
    f(2) = A1*x(0) - (1-x(0))*x(4) - x(2) - theta*x(2)*x(3);

    for ( int_type i = 3; i < n; ++i ) {
      real_type xp2 = 1;
      if ( i+2 < n ) xp2 = x(i+2);
      else if ( i+2 == n ) xp2 = 0;
      f(i) = x(0)*x(i-2) - (1-x(0))*xp2 - x(i) - theta*x(i-1)*x(i);
    }
  }

  virtual
  int_type
  jacobianNnz() const override {
    int_type kk = 10;
    for ( int_type i = 3; i < n; ++i ) {
      kk += 4;
      if ( i+2 < n ) ++kk;
    }
    return kk;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) { ii(kk) = I; jj(kk) = J; ++kk; }

    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,2);

    SETIJ(1,0);
    SETIJ(1,1);
    SETIJ(1,3);

    SETIJ(2,0);
    SETIJ(2,2);
    SETIJ(2,3);
    SETIJ(2,4);

    for ( int_type i = 3; i < n; ++i ) {
      SETIJ(i,i-2);
      SETIJ(i,i-1);
      SETIJ(i,i);
      if ( i+2 < n ) SETIJ(i,i+2);
      SETIJ(i,0);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {

    jac(0) = A0 + x(2);
    jac(1) = - theta*A1;
    jac(2) = x(0) - 1;

    jac(3) = B0 + x(3);
    jac(4) = - theta*A1;
    jac(5) = x(0)-1;

    jac(6) = A1 + x(4);
    jac(7) = -1 - theta*x(3);
    jac(8) =    - theta*x(2);
    jac(9) = x(0)-1;

    int_type kk = 10;

    for ( int_type i = 3; i < n; ++i ) {
      real_type xp2 = 1;
      if      ( i+2 <  n ) xp2 = x(i+2);
      else if ( i+2 == n ) xp2 = 0;
      jac(kk++) = x(0);
      jac(kk++) = -theta*x(i);
      jac(kk++)   = -1 - theta*x(i-1);
      if ( i+2 < n ) jac(kk++) = x(0)-1;
      jac(kk++) = x(i-2)+xp2;
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override
  { }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type istart ) const override {
    for ( int_type i = 0; i < n; ++i ) {
      switch ( i % 8 ) {
        case 0: x(i) = 0.1; break;
        case 1: x(i) = 0.2; break;
        case 2: x(i) = 0.3; break;
        case 3: x(i) = 0.4; break;
        case 4: x(i) = 0.5; break;
        case 5: x(i) = 0.4; break;
        case 6: x(i) = 0.3; break;
        case 7: x(i) = 0.2; break;
      }
    }
    switch ( istart ) {
    case 1: for ( int_type i = 0; i < n; ++i ) x(i) *= 10;  break;
    case 2: for ( int_type i = 0; i < n; ++i ) x(i) *= 100; break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 3; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( std::abs(x(i)) < 10000, "x range" );
  }

};

