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

class Weibull : public nonlinearSystem {
  real_type z[99], y[99];
  int_type const NPT;

public:

  Weibull()
  : nonlinearSystem(
      "Weibull function",
      "@Article{Shan70,\n"
      "  Title   = {Conditioning of quasi-{N}ewton methods for function minimization},\n"
      "  Author  = {David F. Shanno},\n"
      "  Journal = {Mathematics of Computation},\n"
      "  Year    = {1970},\n"
      "  Number  = {111},\n"
      "  Pages   = {647--656},\n"
      "  Volume  = {24}\n"
      "}\n",
      3
    )
  , NPT(99)
  {
    for ( int_type i = 0; i < NPT; ++i ) {
      z[i] = (i+1)*0.01;
      y[i] = 25 + pow( 50 * log(1/z[i]), 2./3. );
    }
  }

  void
  map( dvec_t const & x, dvec_t & eq ) const {
    for ( int_type k = 0; k < NPT; ++k ) {
      real_type y2 = y[k]-x(2);
      real_type g  = pow( y2, x(1) )/ x(0);
      real_type f  = exp( -g );
      eq(k) = f - z[k];
    }
  }

  void
  Grad_map( dvec_t const & x, int_type k, dvec_t & G ) const {
    real_type y2 = y[k]-x(2);
    real_type g  = pow( y2, x(1) )/ x(0);
    real_type fg = g*exp( -g );
    real_type lg = log(y2);
    G(0) = fg/x(0);
    G(1) = -fg*lg;
    G(2) = fg*x(1)/y2;
  }

  void
  Hess_map( dvec_t const & x, int_type k, dmat_t & H ) const {
    real_type y2   = y[k]-x(2);
    real_type g    = pow( y2, x(1) )/ x(0);
    real_type fg   = g*exp( -g );
    real_type fg_1 = (1-g)*exp( -g );
    real_type lg   = log(y2);

    real_type g_x0 = -g/x(0);
    real_type g_x1 = g*lg;
    real_type g_x2 = -g*x(1)/y2;

    //real_type lg_x2 = -1/y2;

    H(0,0) = (fg_1*g_x0-fg/x(0))/x(0);
    H(0,1) = fg_1*g_x1/x(0);
    H(0,2) = fg_1*g_x2/x(0);;

    H(1,0) = H(0,1);
    H(1,1) = -fg_1*g_x1*lg;
    H(1,2) = fg/y2 - fg_1*g_x2*lg;

    H(2,0) = H(0,2);
    H(2,1) = H(1,2);
    H(2,2) = (fg_1*g_x2+fg/y2)*x(1)/y2;
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    dvec_t eq(NPT), G(n);
    map( x, eq );
    f.setZero();
    for ( int_type k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      f += eq(k)*G;
    }
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

  // sum (eq(k,x)-yk) * D_x(j) eq(k,x)
  // sum (eq(k,x)-yk)^2 - z

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    dvec_t eq(NPT), G(n);
    dmat_t H(n,n);
    map( x, eq );
    jac.setZero();
    int_type kk = 0;
    for ( int_type k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      kk = 0;
      for ( int_type i = 0; i < n; ++i ) {
        for ( int_type j = 0; j < n; ++j ) {
          jac(kk) += eq(k)*H(i,j)+G(i)*G(j); ++kk;
        }
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
    x(0) = 250;
    x(1) = 0.3;
    x(2) = 5;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

