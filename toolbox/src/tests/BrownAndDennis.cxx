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

class BrownAndDennis : public nonlinearSystem {
public:

  BrownAndDennis()
  : nonlinearSystem(
    "Brown and Dennis Function",
    "@book{brent2013,\n"
    "  author    = {Brent, R.P.},\n"
    "  title     = {Algorithms for Minimization Without Derivatives},\n"
    "  isbn      = {9780486143682},\n"
    "  series    = {Dover Books on Mathematics},\n"
    "  year      = {2013},\n"
    "  publisher = {Dover Publications}\n"
    "}\n",
    4
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
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);
    
    f(0) = f(1) = f(2) = f(3) = 0;

    for ( int_type i = 0; i < 20; ++i ) {
      real_type c    = (i+1.0) / 5.0;

      real_type sinc = sin(c);
      real_type f1   = x1 + c * x2 - exp(c);
      real_type f2   = x3 + sinc * x4 - cos(c);

      real_type f11 = f1 * f1;
      real_type f22 = f2 * f2;

      real_type tmp = f11+f22;
      real_type t1  = f1*tmp;
      real_type t2  = f2*tmp;

      f(0) += t1;
      f(1) += t1 * c;
      f(2) += t2;
      f(3) += t2 * sinc;
    }

    f(0) *= 4;
    f(1) *= 4;
    f(2) *= 4;
    f(3) *= 4;

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
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    real_type x4 = x(3);

    jac.setZero();

    for ( int_type i = 0; i < 20; ++i ) {

      real_type c = (i+1.0) / 5.0;
      real_type sinc = sin(c);

      real_type f1 = x1 + c * x2 - exp(c);
      real_type f2 = x3 + sinc * x4 - cos(c);

      real_type f11 = f1*f1;
      real_type f12 = f1*f2;
      real_type f22 = f2*f2;

      // f(0) += f1^3+f1*f2^2;
      real_type t1 = 3*f11 + f22;
      real_type t2 = 2*f12;
      jac(0) += t1;
      jac(1) += t1 * c;
      jac(2) += t2;
      jac(3) += t2 * sinc;

      // f(1) += (f1^3+f1*f2^2) * c;
      t1 *= c;
      t2 *= c;
      jac(4) += t1;
      jac(5) += t1 * c;
      jac(6) += t2;
      jac(7) += t2 * sinc;

      // f(2) += f1^2*f2 + f2^3;
      t1 = 2*f12;
      t2 = f11+3*f22;
      jac(8)  += t1;
      jac(9)  += t1 * c;
      jac(10) += t2;
      jac(11) += t2 * sinc;

      // f(3) += (f1^2*f2 + f2^3) * sinc;
      t1 *= sinc;
      t2 *= sinc;
      jac(12) += t1;
      jac(13) += t1 * c;
      jac(14) += t2;
      jac(15) += t2 * sinc;
    }
    jac *= 4;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    x(0) = -11.5844;
    x(1) =  13.1999;
    x(2) = -0.406200;
    x(3) =  0.240998;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 25.0;
    x(1) = 5.0;
    x(2) = -5.0;
    x(3) = -1.0;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  { }

};
