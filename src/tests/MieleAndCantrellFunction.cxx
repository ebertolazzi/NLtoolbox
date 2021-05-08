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

class MieleAndCantrellFunction : public nonlinearSystem {
public:

  MieleAndCantrellFunction()
  : nonlinearSystem(
      "Miele and Cantrell function",
      "@article{Grippo:1991,\n"
      "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
      "  title   = {A Class of Nonmonotone Stabilization Methods\n"
      "             in Unconstrained Optimization},\n"
      "  journal = {Numer. Math.},\n"
      "  year    = {1991},\n"
      "  volume  = {59},\n"
      "  number  = {1},\n"
      "  pages   = {779--805},\n"
      "  doi     = {10.1007/BF01385810},\n"
      "}\n",
      4
    )
  {}

  void
  add_grad1( dvec_t const & x, dvec_t & g ) const {
    real_type tmp = 4*power3( exp(x(0)) - x(1) );
    g(0) += tmp*exp(x(0));
    g(1) -= tmp;
  }

  void
  add_hess1( dvec_t const & x, dmat_t & h ) const {
    real_type e   = exp(x(0));
    real_type t   = e - x(1);
    real_type tte = t*t*e;
    h(0,0) += 4*tte*(4*e-1);
    h(0,1) += -12*tte;
    h(1,0) += -12*tte;
    h(1,1) += 12*t*t;
  }

  void
  add_grad2( dvec_t const & x, dvec_t & g ) const {
    real_type tmp = 600*power5(x(1)-x(2));
    g(1) += tmp;
    g(2) -= tmp;
  }

  void
  add_hess2( dvec_t const & x, dmat_t & h ) const {
    real_type tmp = 3000*power4(x(1)-x(2));
    h(1,1) += tmp;
    h(1,2) -= tmp;
    h(2,1) -= tmp;
    h(2,2) += tmp;
  }

  void
  add_grad3( dvec_t const & x, dvec_t & g ) const {
    real_type t   = tan(x(2)-x(3));
    real_type tmp = 4*power3(t)*(1+power2(t));
    g(2) += tmp;
    g(3) -= tmp;
  }

  void
  add_hess3( dvec_t const & x, dmat_t & h ) const {
    real_type t   = tan(x(2)-x(3));
    real_type t2  = t*t;
    real_type tmp = ((20*t2+32)*t2+12)*t2;
    h(2,2) += tmp;
    h(2,3) -= tmp;
    h(3,2) -= tmp;
    h(3,3) += tmp;
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
    f(0) = 8*x(0)*power6(x(0));
    f(1) = f(2) = f(3) = 0;
    add_grad2( x, f );
    add_grad3( x, f );
    add_grad1( x, f );
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
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    dmat_t h(4,4);
    h.setZero();
    h(0,0) = 8*7*power6(x(0));
    add_hess2( x, h );
    add_hess3( x, h );
    add_hess1( x, h );
    int_type kk = 0;
    for ( int_type i = 0; i < 4; ++i )
      for ( int_type j = 0; j < 4; ++j )
        jac(kk++) = h(i,j);
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 1;
    x(2) = 1;
    x(3) = 1;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type init_point ) const override {
    switch ( init_point ) {
    case 0:
      x(0) = 10;
      x(1) = -10;
      x(2) = -10;
      x(3) = -10;
    break;
    case 1:
      x(0) = 1;
      x(1) = 2;
      x(2) = 2;
      x(3) = 2;
    break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 2; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
