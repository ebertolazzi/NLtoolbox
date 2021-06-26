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

class ChebyquadFunction : public nonlinearSystem {

  mutable real_type T[10];
  mutable real_type dT[10];

public:
  // Only the values N = 1, 2, 3, 4, 5, 6, 7 and 9 may be used.
  ChebyquadFunction( integer dim )
  : nonlinearSystem(
      "Chebyquad function",
      "@article {Fletcher:1965,\n"
      "  author  = {Fletcher, R.},\n"
      "  title   = {Function minimization without evaluating derivatives -- {A} review},\n"
      "  journal = {The Computer Journal},\n"
      "  year    = {1965},\n"
      "  volume  = {8},\n"
      "  pages   = {33--41},\n"
      "  doi     = {10.1093/comjnl/8.1.33},\n"
      "}\n\n"
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
      dim
    )
  {
    UTILS_ASSERT(
      dim > 0 && dim < 10,
      "ChebyquadFunction:: dimension n = {} must be on [1,2,...,8]", dim
    );
  }

  void
  Chebyshev( real_type x ) const {
    T[0] = 1;
    T[1] = 2*x-1;
    for ( integer j=1; j < n; ++j )
      T[j+1] = 2*x*T[j]-T[j-1];
  }

  void
  Chebyshev_D( real_type x ) const {
    T[0]  = 1;
    T[1]  = 2*x-1;
    dT[0] = 0;
    dT[1] = 2;
    for ( integer j=1; j < n; ++j ) {
      T[j+1]  = 2*x*T[j] - T[j-1];
      dT[j+1] = 2*x*dT[j] - dT[j-1] + 2*T[j];
    }
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    real_type f = 0;
    for ( integer j = 0; j < n; ++j ) {
      Chebyshev( x(j) );
      f += T[k+1];
    }
    f /= real_type(n);
    if ( (k%2) == 1 ) f += 1.0/(power2(k+1)-1);
    return f;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer k = 0; k < n; ++k ) f(k) = 0;
    for ( integer j = 0; j < n; ++j ) {
      Chebyshev( x(j) );
      for ( integer i = 0; i < n; ++i )
        f(i) += T[i+1];
    }
    for ( integer k = 0; k < n; ++k ) {
      f(k) /= real_type(n);
      if ( (k%2) == 1 ) f(k) += 1.0/(power2(k+1)-1);
    }
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer j = 0; j < n; ++j )
      for ( integer i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer j = 0; j < n; ++j ) {
      Chebyshev_D( x(j) );
      for ( integer i = 0; i < n; ++i )
        jac(kk++) = dT[i+1]/n;
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    switch (n) {
    case 1:
      x << 0.5;
      break;
    case 2:
      x << 0.2113248654051871, 0.7886751345948129;
      break;
    case 3:
      x << 0.1464466094067262, 0.5, 0.8535533905932737;
      break;
    case 4:
      x << 0.1026727638541169, 0.40620376295746,
           0.5937962370425399, 0.8973272361458831;
      break;
    case 5:
      x << 0.08375125649950906, 0.3127292952232095, 0.5,
           0.6872707047767905, 0.9162487435004909;
      break;
    case 6:
      x << 0.06687659094608972, 0.3666822992416476, 0.2887406731194442,
           0.7112593268805557,  0.6333177007583524,  0.9331234090539103;
      break;
    case 7:
      x << 0.0580691496209755, 0.2351716123574216, 0.3380440947400462, 0.5,
           0.6619559052599538, 0.7648283876425784, 0.9419308503790246;
      break;
    case 9:
      x << 0.04420534613578277, 0.199490672309881,  0.23561910847106,
           0.4160469078925981,  0.4999999999999999, 0.5839530921074021,
           0.76438089152894,    0.8005093276901191, 0.9557946538642172;
      break;
    }
  }

  integer
  numExactSolution() const override {
    if ( n < 10 && n > 0 && n != 8 ) return 1;
    return 0;
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    for ( integer i = 0; i < n; ++i ) {
      real_type s = (i+1)/real_type(n+1);
      x(i) = s;
    }
  }

  integer
  numInitialPoint() const override { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for ( integer i = 0; i < n; ++i )
    //  UTILS_ASSERT(:abs(x(i)) < 5000, "x range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    //for ( integer i = 0; i < n; ++i )
    //  { U[i] = 5000; L[i] = -5000; }
    U.fill(real_max);
    L.fill(-real_max);
  }

};
