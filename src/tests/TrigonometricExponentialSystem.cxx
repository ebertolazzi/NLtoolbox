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

#define TRIGONOMETRIC_EXPONENTIAL_BIBTEX \
"Spedicato, E.\n" \
"Computational experience with quasi-newton algoritms.\n" \
"for minimization problems of moderately large size.\n" \
"Rep. CISE-N-175, Segrate (Milano), 1975.\n\n" \
"@article{More:1981,\n" \
"  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n" \
"  title   = {Testing Unconstrained Optimization Software},\n" \
"  journal = {ACM Trans. Math. Softw.},\n" \
"  year    = {1981},\n" \
"  volume  = {7},\n" \
"  number  = {1},\n" \
"  pages   = {17--41},\n" \
"  doi     = {10.1145/355934.355936},\n" \
"}\n"

class TrigonometricExponentialSystem1 : public nonlinearSystem {
public:
  
  TrigonometricExponentialSystem1( integer neq )
  : nonlinearSystem(
      "Trigonometric Exponential System prob 1",
      TRIGONOMETRIC_EXPONENTIAL_BIBTEX,
      neq)
  { checkEven(n,2); }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( (i%2) == 0 ) return 3*power3(x(i))+2*x(i+1)-5 + sin(x(i)-x(i+1))*sin(x(i)+x(i+1));
    return -x(i-1)*exp(x(i-1)-x(i)) + 4*x(i)-3;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n; i+=2 )
      f(i) = 3*power3(x(i))+2*x(i+1)-5 + sin(x(i)-x(i+1))*sin(x(i)+x(i+1));
    for ( integer i = 1; i < n; i+=2 )
      f(i) = -x(i-1)*exp(x(i-1)-x(i)) + 4*x(i)-3;
  }

  integer
  jacobianNnz() const override {
    return 2*n;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( integer i = 0; i < n; i+=2 ) {
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    for ( integer i = 1; i < n; i+=2 ) {
      SETIJ(i,i);
      SETIJ(i,i-1);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; i+=2 ) {
      jac(kk++) = 9*power2(x(i)) + sin(2*x(i));
      jac(kk++) = 2 - sin(2*x(i+1));
    }
    for ( integer i = 1; i < n; i+=2 ) {
      jac(kk++) = x(i-1)*exp(x(i-1)-x(i)) + 4;
      jac(kk++) = -(x(i-1)+1)*exp(x(i-1)-x(i));
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }
  
  void
  getInitialPoint( dvec_t & x, integer ) const override {
    for ( integer k = 0; k < n; ++k ) x(k) = 0;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT0( abs(x(i)) < 100, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(100);
    L.fill(-100);
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class TrigonometricExponentialSystem2 : public nonlinearSystem {
public:
  
  TrigonometricExponentialSystem2( integer neq )
  : nonlinearSystem(
      "Trigonometric Exponential System prob 2",
      TRIGONOMETRIC_EXPONENTIAL_BIBTEX,
      neq
    )
  { checkOdd(n,6); }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 3*power3(x(0)-x(2))
         + 2*x(1)-5
         + sin( x(0)-x(1)-x(2) )*sin( x(0)+x(1)-x(2) );

    for ( integer i = 1; i < n-1; i += 2 )
      f(i) = (x(i+1)-x(i-1))*exp(x(i-1)-x(i)-x(i+1))+4*x(i) - 3;

    for ( integer i = 2; i < n-1; i += 2 )
      f(i) = 3*power3(x(i)-x(i+2)) + 6*power3(x(i)-x(i-2)) +
             2*x(i+1)-4*x(i-1)+5
             -2*sin( x(i-2)-x(i-1)-x(i) )*sin( x(i-2)+x(i-1)-x(i) )
             +sin( x(i)-x(i+1)-x(i+2) )*sin( x(i)+x(i+1)-x(i+2) );

    f(n-1) = - 6*power3(x(n-1)-x(n-3))
             - 4*x(n-2) + 10
             - 2*sin( x(n-3)-x(n-2)-x(n-1) )*sin( x(n-3)+x(n-2)-x(n-1) );
  }

  integer
  jacobianNnz() const override {
    integer kk = 6;
    for ( integer i = 1; i < n-1; i += 2 ) kk += 3;
    for ( integer i = 2; i < n-1; i += 2 ) kk += 5;
    return kk;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,2);

    for ( integer i = 1; i < n-1; i += 2 ) {
      SETIJ(i,i-1);
      SETIJ(i,i);
      SETIJ(i,i+1);
    }

    for ( integer i = 2; i < n-1; i += 2 ) {
      SETIJ(i,i-2);
      SETIJ(i,i-1);
      SETIJ(i,i-0);
      SETIJ(i,i+1);
      SETIJ(i,i+2);
    }

    SETIJ(n-1,n-3);
    SETIJ(n-1,n-2);
    SETIJ(n-1,n-1);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    jac(kk++) = 9*power2(x(0)-x(2))+sin(2*(x(0)-x(2)));
    jac(kk++) = 2-sin(2*x(1));
    jac(kk++) = -9*power2(x(0)-x(2))-sin(2*(x(0)-x(2)));

    for ( integer i = 1; i < n-1; i += 2 ) {
      real_type ex = exp(x(i-1)-x(i)-x(i+1));
      real_type tp = x(i+1)-x(i-1)-1;
      jac(kk++) = tp*ex;
      jac(kk++) = 4+(x(i-1)-x(i+1))*ex;
      jac(kk++) = -tp*ex;
    }

    for ( integer i = 2; i < n-1; i += 2 ) {
      jac(kk++) = 2*sin(2*(x(i)-x(i-2)))-18*power2(x(i)-x(i-2));
      jac(kk++) = -4+2*sin(2*x(i-1));
      jac(kk++) = sin(2*(x(i)-x(i+2)))-2*sin(2*(x(i)-x(i-2)))
                + 27*power2(x(i))
                - 36*x(i)*x(i-2)
                - 18*x(i)*x(i+2)
                + 18*power2(x(i-2))
                + 9*power2(x(i+2));
      jac(kk++) = 2-sin(2*x(i+1));
      jac(kk++) = -9*power2(x(i)-x(i+2))-sin(2*(x(i)-x(i+2)));
    }

    jac(kk++) = -2*sin(2*(x(n-3)-x(n-1)))+18*power2(x(n-3)-x(n-1));
    jac(kk++) = 2*sin(2*x(n-2))-4;
    jac(kk++) = 2*sin(2*(x(n-3)-x(n-1)))-18*power2(x(n-3)-x(n-1));
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }
  
  void
  getInitialPoint( dvec_t & x, integer ) const override {
    for ( integer k = 0; k < n; ++k ) x(k) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT0( abs(x(i)) < 100, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(100);
    L.fill(-100);
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/


class TrigExp : public nonlinearSystem {
public:
  
  TrigExp( integer neq )
  : nonlinearSystem(
      "TrigExp",
      "@article{Ruggiero:1992,\n"
      "  author  = {Gomes-Ruggiero, M. and Martínez, J. and Moretti, A.},\n"
      "  title   = {Comparing Algorithms for Solving Sparse Nonlinear\n"
      "             Systems of Equations},\n"
      "  journal = {SIAM Journal on Scientific and Statistical Computing},\n"
      "  volume  = {13},\n"
      "  number  = {2},\n"
      "  pages   = {459-483},\n"
      "  year    = {1992},\n"
      "  doi     = {10.1137/0913025},\n"
      "}\n",
      neq
    )
  { checkMinEquations(n,3); }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( i == 0 )
      return 3*x(0)*x(0)+2*x(1)-5 + sin(x(0)-x(1))*sin(x(0)+x(1));
    else if ( i == n-1 )
      return -x(n-2)*exp(x(n-1)-x(n-2)) + 4*x(n-1)-3;
    return -x(i-1)*exp(x(i-1)-x(i))
           + x(i)*(4+3*x(i)*x(i))+2*x(i+1)
           + sin(x(i)-x(i+1))*sin(x(i)+x(i+1));
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0)   = 3*x(0)*x(0)+2*x(1)-5 + sin(x(0)-x(1))*sin(x(0)+x(1));
    f(n-1) = -x(n-2)*exp(x(n-1)-x(n-2)) + 4*x(n-1)-3;
    for ( integer i = 1; i < n-1; ++i )
      f(i) = -x(i-1)*exp(x(i-1)-x(i))
             + x(i)*(4+3*x(i)*x(i))+2*x(i+1)
             + sin(x(i)-x(i+1))*sin(x(i)+x(i+1));
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
    SETIJ(0,1);
    SETIJ(n-1,n-1);
    SETIJ(n-1,n-2);
    for ( integer i = 1; i < n-1; ++i ) {
      SETIJ(i,i-1);
      SETIJ(i,i);
      SETIJ(i,i+1);
    }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    jac(kk++) = 6*x(0) + cos(x(0)-x(1))*sin(x(0)+x(1))
              + sin(x(0)-x(1))*cos(x(0)+x(1));
    jac(kk++) = 2 - cos(x(0)-x(1))*sin(x(0)+x(1))
              + sin(x(0)-x(1))*cos(x(0)+x(1));
    jac(kk++) = 4 - x(n-2)*exp(x(n-1)-x(n-2));
    jac(kk++) = (x(n-2)-1)*exp(x(n-1)-x(n-2));
    for ( integer i = 1; i < n-1; ++i ) {
      jac(kk++) = -(1+x(i-1))*exp(x(i-1)-x(i));
      jac(kk++) = x(i-1)*exp(x(i-1)-x(i))+9*x(i)*x(i)+4
                + cos(x(i)-x(i+1))*sin(x(i)+x(i+1))
                + sin(x(i)-x(i+1))*cos(x(i)+x(i+1));
      jac(kk++) = 2-cos(x(i)-x(i+1))*sin(x(i)+x(i+1))
                + sin(x(i)-x(i+1))*cos(x(i)+x(i+1));
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }
  
  void
  getInitialPoint( dvec_t & x, integer ini ) const override {
    switch( ini ) {
    case 0:
      for ( integer k = 0; k < n; ++k ) x(k) = 1 / real_type(n);
      break;
    case 1:
      for ( integer k = 0; k < n; ++k ) x(k) = 0;
      break;
    case 2:
      for ( integer k = 0; k < n; ++k ) x(k) = 0.3;
      break;
    }
  }

  integer
  numInitialPoint() const override
  { return 3; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT0( abs(x(i)) < 1000, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(1000);
    L.fill(-1000);
  }

};
