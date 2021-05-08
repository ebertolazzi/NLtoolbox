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

class Gulf : public nonlinearSystem {
  real_type rr[99];
public:

  Gulf()
  : nonlinearSystem(
      "Gulf R&D Function",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      3
    )
  {
    for ( int_type i = 0; i < 99; ++i ) {
      real_type arg = (i+1) / 100.0;
      rr[i] = pow(-50.0*log(arg),2.0/3.0) + 25.0;
    }
  }

  real_type
  t_fun( dvec_t const & x, int_type i ) const {
    real_type arg = (i+1) / 100.0;
    real_type   r = rr[i] - x(1);
    if ( r <= 0 ) return nan("t_fun");
    return exp(-pow(r,x(2))/x(0)) - arg;
  }

  void
  t_grad( dvec_t const & x, int_type i, dvec_t & g ) const {
    real_type t1 = rr[i]-x(1);
    if ( t1 <= 0 ) {
      g(0) = g(1) = g(2) = nan("t_grad");
      return;
    }
    real_type t2  = pow(t1,x(2));
    real_type t3  = x(0)*x(0);
    real_type t6  = 1/x(0);
    real_type t8  = exp(-t6*t2);
    real_type t15 = log(t1);
    g(0) =  t8/t3*t2;
    g(1) =  t8*t6/t1*x(2)*t2;
    g(2) = -t8*t6*t15*t2;
  }

  void
  t_hess( dvec_t const & x, int_type i, dmat_t & h ) const {
    real_type t1 = rr[i]-x(1);
    if ( t1 <= 0 ) {
      h.fill(nan("t_hess"));
      return;
    }
    real_type t2  = pow(t1,x(2));
    real_type t5  = exp(-t2/x(0));
    real_type t6  = t5*t2;
    real_type t9  = x(0)*x(0);
    real_type t10 = t9*t9;
    real_type t15 = t5*x(2)*t2;
    real_type t16 = t2-x(0);
    real_type t18 = 1/t9/x(0);
    real_type t20 = 1/t1;
    real_type t23 = log(t1);
    real_type t30 = t1*t1;
    real_type t33 = 1/t9;
    real_type t42 = t23*t23;
    h(0,0) = (t2-2.0*x(0))*t6/t10;
    h(0,1) = t20*t18*t16*t15;
    h(0,2) = -t18*t16*t5*t23*t2;
    h(1,1) = t33/t30*(x(2)*t16+x(0))*t15;
    h(1,2) = -t33*t20*(x(2)*t23*t16-x(0))*t6;
    h(2,2) = t33*t16*t5*t42*t2;
    h(1,0) = h(0,1);
    h(2,0) = h(0,2);
    h(2,1) = h(1,2);
  }

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
    dvec_t g(3);
    f(0) = f(1) = f(2) = 0;
    for ( int_type i = 0; i < 99; ++i ) {
      real_type t = t_fun( x, i );
      t_grad( x, i, g );
      f(0) += t*g(0);
      f(1) += t*g(1);
      f(2) += t*g(2);
    }
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 9;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,2);
    SETIJ(1,1);
    SETIJ(1,2);
    SETIJ(2,2);
    SETIJ(1,0);
    SETIJ(2,0);
    SETIJ(2,1);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = jac(1) = jac(2) = jac(3) = jac(4) = jac(5) = 0;
    for ( int_type i = 0; i < 99; ++i ) {
      dvec_t g(3);
      dmat_t h(3,3);
      real_type t = t_fun( x, i );
      t_grad( x, i, g );
      t_hess( x, i, h );
      jac(0) += t*h(0,0)+g(0)*g(0);
      jac(1) += t*h(0,1)+g(0)*g(1);
      jac(2) += t*h(0,2)+g(0)*g(2);

      jac(3) += t*h(1,1)+g(1)*g(1);
      jac(4) += t*h(1,2)+g(1)*g(2);

      jac(5) += t*h(2,2)+g(2)*g(2);
    }
    jac(6) = jac(1);
    jac(7) = jac(2);
    jac(8) = jac(4);
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 50.0;
    x(1) = 25.0;
    x(2) = 1.5;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 40.0;
    x(1) = 20.0;
    x(2) = 1.20;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    NONLIN_ASSERT( x(0) > 0, "x0!!!!" );
    NONLIN_ASSERT( x(1) > 0, "x0!!!!" );
    for ( int_type i = 0; i < 99; ++i ) {
      real_type t1 = rr[i]-x(1);
      NONLIN_ASSERT( t1 >= 0, "r < 0!!!! r = " << t1 );
    }
  }

  virtual
  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(real_max);
    L.fill(-real_max);
    L[0] = L[1] = 0;
  }

};
