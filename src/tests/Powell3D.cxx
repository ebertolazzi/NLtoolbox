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

class Powell3D : public nonlinearSystem {
public:

  Powell3D()
  : nonlinearSystem(
      "Powell 3D Function",
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
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  virtual
  void
  evalF( dvec_t const & X, dvec_t & f ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);

    real_type t3 = x-y;
    real_type t6 = t3*t3;
    real_type t8 = pow(1.0+t6,2.0);
    real_type t13 = 2.0*t3/t8;
    real_type t17 = cos(m_pi*y*z/2.0);
    real_type t18 = m_pi*t17;
    f(0) = t13;
    f(1) = -t13-z*t18/2.0;
    f(2) = -y*t18/2.0;

    if ( y != 0 ) {
      real_type tt = exp(-(x+2.0*y+z)/y)/y;
      f(0) -= tt;
      f(1) += (x+z)*tt/y;
      f(2) -= tt;
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
    SETIJ(1,0);
    SETIJ(1,1);
    SETIJ(1,2);
    SETIJ(2,0);
    SETIJ(2,1);
    SETIJ(2,2);

    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);

    real_type t3 = x*x;
    real_type t7 = x*y;
    real_type t9 = y*y;
    real_type t13 = -6.0*t3+12.0*t7-6.0*t9+2.0;
    real_type t15 = t3-2.0*t7+t9+1.0;
    real_type t16 = t15*t15;
    real_type t18 = 1/t16/t15;
    real_type t20 = -t18*t13;
    real_type t22 = pow(x-y,2.0);
    real_type t26 = pow(1.0+t22,2.0);
    real_type t29 = m_pi*m_pi;
    real_type t30 = z*z;
    real_type t34 = m_pi*y*z/2.0;
    real_type t35 = sin(t34);
    real_type t42 = cos(t34);
    real_type t46 = (y*t35*m_pi*z-2.0*t42)*m_pi/4.0;

    jac(0) = t18*t13;
    jac(1) = t20;
    jac(2) = 0.0;
    jac(3) = t20;
    jac(4) = -8.0*t18*t22+2.0/t26+t35*t30*t29/4.0;
    jac(5) = t46;
    jac(6) = 0.0;
    jac(7) = t46;
    jac(8) = t35*t9*t29/4.0;

    if ( y != 0 ) {
      real_type t2  = 2.0*y;
      real_type t8  = exp((-x-t2-z)/y);
      real_type t10 = y*y;
      real_type t12 = t8/t10;
      real_type t17 = t8*(-y+x+z)/t10/y;
      real_type t21 = t10*t10;
      jac(0) += t12;
      jac(1) += -t17;
      jac(2) += t12;
      jac(3) += -t17;
      jac(4) += (x+z)*t8*(x+z-t2)/t21;
      jac(5) += -t17;
      jac(6) += t12;
      jac(7) += -t17;
      jac(8) += t12;
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 1;
    x(2) = 1;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 1;
    x(2) = 2;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & X ) const override {
    //real_type x = x(0);
    //real_type y = x(1);
    //real_type z = x(2);
    //ASSERT( y > 0, "X!!!!" );
    //for (  i = 0; i < n; ++i ) {
    //  NONLIN_ASSERT( std::abs(x(i)) < 100, "X!!!!" );
    //}
  }

};
