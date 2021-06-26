/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*
  (1/2) sum (eq(k,x)-yk)^2

  -->

  sum (eq(k,x)-yk) * D_x(j) eq(k,x)
 
*/

/*
  (1/2) sum (e(k)-yk)^2 + sum l(k)*(e(k)-eq(k,x))

  e(k)-zk+l(k) = 0
  e(k)-eq(k,x) = 0
  -sum_k D_x(j) eq(k,x) * l(k) = 0

  // si può eliminare l(k)

  eq(k,x)-e(k) = 0
  sum_k D_x(j) eq(k,x) * (e(k)-yk) = 0    j = 1,2,.,,nx
*/

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

#define BIGGS_BIBTEX \
"@Article{Biggs:1971,\n" \
"  author = {M. C. Biggs},\n" \
"  title  = { Minimization algorithms making use of non-quadratic\n" \
"             properties of the objective function},\n" \
"  volume = {8},\n" \
"  pages  = {315--327},\n" \
"  year   = {1971},\n" \
"  journal = {Journal of the Institute of Mathematics and its Applications}\n" \
"}\n\n" \
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

class BiggsEXP2function : public nonlinearSystem {
  dvec_t z, y;
  integer const NPT;

public:

  BiggsEXP2function()
  : nonlinearSystem( "Biggs EXP2 function", BIGGS_BIBTEX, 2 )
  , NPT(10)
  {
    z.resize(NPT);
    y.resize(NPT);
    for ( integer i = 0; i < NPT; ++i ) z(i) = (i+1)*0.1;
    y = exp(-z.array())-5*exp(-10*z.array());
  }

  void
  map( dvec_t const & x, dvec_t & eq ) const {
    eq = exp(-x(0)*z.array()) - 5*exp(-x(1)*z.array()) - y.array();
  }

  void
  Grad_map( dvec_t const & x, integer k, dvec_t & G ) const {
    G(0) = -z(k)*exp(-x(0)*z(k));
    G(1) = 5*z(k)*exp(-x(1)*z(k));
  }

  void
  Hess_map( dvec_t const & x, integer k, dmat_t & H ) const {
    real_type zk  = z(k);
    real_type zk2 = zk*zk;
    real_type ex0 = exp(-x(0)*zk);
    real_type ex1 = exp(-x(1)*zk);
    H(0,0) = zk2*ex0;
    H(0,1) = 0;
    H(1,1) = -5*zk2*ex1;
    H(1,0) = 0;
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    dvec_t eq(NPT), G(n);
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      f += eq(k)*G;
    }
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j )
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
    integer kk = 0;
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      kk = 0;
      for ( integer i = 0; i < n; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          jac(kk) += eq(k)*H(i,j)+G(i)*G(j); ++kk;
        }
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = 16.7046761257;
    x(1) = 16.7046761257;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 2;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP3function : public nonlinearSystem {
  dvec_t z, y;
  integer const NPT;

public:

  BiggsEXP3function()
  : nonlinearSystem( "Biggs EXP3 function", BIGGS_BIBTEX, 3 )
  , NPT(10)
  {
    z.resize(NPT);
    y.resize(NPT);
    for ( integer i = 0; i < NPT; ++i ) z(i) = (i+1)*0.1;
    y = exp(-z.array())-5*exp(-10*z.array());
  }

  void
  map( dvec_t const & x, dvec_t & eq ) const {
    eq = exp(-x(0)*z.array()) - x(2)*exp(-x(1)*z.array()) - y.array();
  }

  void
  Grad_map( dvec_t const & x, integer k, dvec_t & G ) const {
    real_type zk  = z(k);
    real_type ex0 = exp(-x(0)*zk);
    real_type ex1 = exp(-x(1)*zk);
    G(0) = -zk*ex0;
    G(1) = x(2)*zk*ex1;
    G(2) = -ex1;
  }

  void
  Hess_map( dvec_t const & x, integer k, dmat_t & H ) const {
    real_type zk  = z(k);
    real_type zk2 = zk*zk;
    real_type ex0 = exp(-x(0)*zk);
    real_type ex1 = exp(-x(1)*zk);
    H(0,0) = zk2*ex0;  H(0,1) = 0;             H(0,2) = 0;
    H(1,0) = 0;        H(1,1) = -x(2)*zk2*ex1; H(1,2) = zk*ex1;
    H(2,0) = 0;        H(2,1) = zk*ex1;        H(2,2) = 0;
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    dvec_t eq(NPT), G(n);
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      f += eq(k)*G;
    }
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j )
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
    integer kk = 0;
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      kk = 0;
      for ( integer i = 0; i < n; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          jac(kk) += eq(k)*H(i,j)+G(i)*G(j); ++kk;
        }
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = 1;
    x(1) = 1;
    x(2) = 5;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 2;
    x(2) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  k = 0; k < n; ++k )
    //  UTILS_ASSERT( abs(x(k)) < 1000, "Bad Range" );
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP4function : public nonlinearSystem {
  dvec_t z, y;
  integer const NPT;

public:

  BiggsEXP4function()
  : nonlinearSystem( "Biggs EXP4 function", BIGGS_BIBTEX, 4 )
  , NPT(10)
  {
    z.resize(NPT);
    y.resize(NPT);
    for ( integer i = 0; i < NPT; ++i ) z(i) = (i+1)*0.1;
    y = exp(-z.array())-5*exp(-10*z.array());
  }

  void
  map( dvec_t const & x, dvec_t & eq ) const {
    eq = x(2)*exp(-x(0)*z.array()) - x(3)*exp(-x(1)*z.array()) - y.array();
  }

  void
  Grad_map( dvec_t const & x, integer k, dvec_t & G ) const {
    real_type zk  = z(k);
    real_type ex0 = exp(-x(0)*zk);
    real_type ex1 = exp(-x(1)*zk);
    G(0) = -x(2)*zk*ex0;
    G(1) = x(3)*zk*ex1;
    G(2) = ex0;
    G(3) = -ex1;
  }

  void
  Hess_map( dvec_t const & x, integer k, dmat_t & H ) const {
    real_type zk  = z(k);
    real_type zk2 = zk*zk;
    real_type ex0 = exp(-x(0)*zk);
    real_type ex1 = exp(-x(1)*zk);

    H(0,0) = x(2)*zk2*ex0;
    H(0,1) = 0;
    H(0,2) = -zk*ex0;
    H(0,3) = 0;

    H(1,0) = 0;
    H(1,1) = -x(3)*zk2*ex1;
    H(1,2) = 0;
    H(1,3) = zk*ex1;

    H(2,0) = -zk*ex0;
    H(2,1) = 0;
    H(2,2) = 0;
    H(2,3) = 0;

    H(3,0) = 0;
    H(3,1) = zk*ex1;
    H(3,2) = 0;
    H(3,3) = 0;
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    dvec_t eq(NPT), G(n);
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      f += eq(k)*G;
    }
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j )
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
    integer kk = 0;
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      kk = 0;
      for ( integer i = 0; i < n; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          jac(kk) += eq(k)*H(i,j)+G(i)*G(j); ++kk;
        }
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = -1.370515321;
    x(1) = -1.370515321;
    x(2) = 0.146554089168054;
    x(3) = 0.003045452799256;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 2;
    x(2) = 1;
    x(3) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  k = 0; k < n; ++k )
    //  UTILS_ASSERT( abs(x(k)) < 1000, "Bad Range" );
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP5function : public nonlinearSystem {
  dvec_t z, y;
  integer const NPT;

public:

  BiggsEXP5function()
  : nonlinearSystem( "Biggs EXP5 function", BIGGS_BIBTEX, 5 )
  , NPT(11)
  {
    z.resize(NPT);
    y.resize(NPT);
    for ( integer k = 0; k < NPT; ++k ) z(k) = (k+1)*0.1;
    y = exp(-z.array())-5*exp(-10*z.array()) + 3*exp(-4*z.array());
  }

  void
  map( dvec_t const & x, dvec_t & eq ) const {
    eq = x(2)*exp(-x(0)*z.array())
       - x(3)*exp(-x(1)*z.array())
       + 3*exp(-x(4)*z.array())
       - y.array();
  }

  void
  Grad_map( dvec_t const & x, integer k, dvec_t & G ) const {
    real_type zk  = z(k);
    real_type ex0 = exp(-x(0)*zk);
    real_type ex1 = exp(-x(1)*zk);
    real_type ex4 = exp(-x(4)*zk);
    G(0) = -x(2)*zk*ex0;
    G(1) = x(3)*zk*ex1;
    G(2) = ex0;
    G(3) = -ex1;
    G(4) = -3*zk*ex4;
  }

  void
  Hess_map( dvec_t const & x, integer k, dmat_t & H ) const {
    real_type zk  = z(k);
    real_type zk2 = zk*zk;
    real_type ex0 = exp(-x(0)*zk);
    real_type ex1 = exp(-x(1)*zk);
    real_type ex4 = exp(-x(4)*zk);

    H(0,0) = x(2)*zk2*ex0;
    H(0,1) = 0;
    H(0,2) = -zk*ex0;
    H(0,3) = 0;
    H(0,4) = 0;

    H(1,0) = 0;
    H(1,1) = -x(3)*zk2*ex1;
    H(1,2) = 0;
    H(1,3) = zk*ex1;
    H(1,4) = 0;

    H(2,0) = -zk*ex0;
    H(2,1) = 0;
    H(2,2) = 0;
    H(2,3) = 0;
    H(2,4) = 0;

    H(3,0) = 0;
    H(3,1) = zk*ex1;
    H(3,2) = 0;
    H(3,3) = 0;
    H(3,4) = 0;

    H(4,0) = 0;
    H(4,1) = 0;
    H(4,2) = 0;
    H(4,3) = 0;
    H(4,4) = 3*zk2*ex4;
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    dvec_t eq(NPT), G(n);
    map( x, eq );
    f.setZero();
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      f += eq(k)*G;
    }
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j )
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
    integer kk = 0;
    for ( integer k = 0; k < NPT; ++k ) {
      Grad_map( x, k, G );
      Hess_map( x, k, H );
      kk = 0;
      for ( integer i = 0; i < n; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          jac(kk) += eq(k)*H(i,j)+G(i)*G(j); ++kk;
        }
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx ) const override {
    x(0) = 2.645310671110160;
    x(1) = 2.693935021162610;
    x(2) = 0.048826924681352;
    x(3) = 0.049142366883430;
    x(4) = 0.025390267672614;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 2;
    x(2) = 1;
    x(3) = 1;
    x(4) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  k = 0; k < n; ++k )
    //  UTILS_ASSERT( abs(x(k)) < 1000, "Bad Range" );
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class BiggsEXP6function : public nonlinearSystem {
  real_type const xmin;
  real_type const xmax;
  dvec_t y, t;
public:

  BiggsEXP6function()
  : nonlinearSystem( "Biggs EXP6 function", BIGGS_BIBTEX, 6 )
  , xmin(-10)
  , xmax(20)
  {
    y.resize(6);
    t.resize(6);
    for ( integer i = 0; i < 6; ++i ) t(i) = (i+1)*0.1;
    y = exp(-t.array())-5*exp(-10*t.array())+3*exp(-4*t.array());
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    real_type e0 = exp(-t(k)*x(0));
    real_type e1 = exp(-t(k)*x(1));
    real_type e4 = exp(-t(k)*x(4));
    return x(2)*e0 -x(3)*e1 +x(5)*e4 - y(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer k = 0; k < 6; ++k ) {
      real_type e0 = exp(-t(k)*x(0));
      real_type e1 = exp(-t(k)*x(1));
      real_type e4 = exp(-t(k)*x(4));
      f(k) = x(2)*e0 - x(3)*e1 + x(5)*e4 - y(k);
    }
  }

  integer
  jacobianNnz() const override
  { return 36; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer k = 0; k < 6; ++k ) {
      ii(kk) = k; jj(kk) = 0; ++kk;
      ii(kk) = k; jj(kk) = 1; ++kk;
      ii(kk) = k; jj(kk) = 2; ++kk;
      ii(kk) = k; jj(kk) = 3; ++kk;
      ii(kk) = k; jj(kk) = 4; ++kk;
      ii(kk) = k; jj(kk) = 5; ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer k = 0; k < 6; ++k ) {
      real_type e0 = exp(-t(k)*x(0));
      real_type e1 = exp(-t(k)*x(1));
      real_type e4 = exp(-t(k)*x(4));
      jac(kk) = -x(2)*t(k)*e0; ++kk;
      jac(kk) = x(3)*t(k)*e1;  ++kk;
      jac(kk) = e0;            ++kk;
      jac(kk) = -e1;           ++kk;
      jac(kk) = -x(5)*t(k)*e4; ++kk;
      jac(kk) =  e4;           ++kk;
    }
  }

  void
  getExactSolution( dvec_t & x, integer idx) const override {
    switch ( idx ) {
    case 0:
      x(0) = 10;
      x(1) = 4;
      x(2) = -5;
      x(3) = -3;
      x(4) = 1;
      x(5) = 1;
      break;
    case 1:
      x(0) = 1;
      x(1) = 10;
      x(2) = 1;
      x(3) = 5;
      x(4) = 4;
      x(5) = 3;
      break;
    }
  }

  integer
  numExactSolution() const override
  { return 2; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 2;
    x(2) = 1;
    x(3) = 1;
    x(4) = 1;
    x(5) = 1;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    UTILS_ASSERT0(
      x(0) > xmin && x(0) < xmax &&
      x(1) > xmin && x(1) < xmax &&
      x(2) > xmin && x(2) < xmax &&
      x(3) > xmin && x(3) < xmax &&
      x(4) > xmin && x(4) < xmax &&
      x(5) > xmin && x(5) < xmax,
      "Bad Range"
    );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    L.fill(xmin);
    U.fill(xmax);
  }

};
