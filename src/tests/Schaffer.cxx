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

#define SHAFFER_BIBTEX \
"@book{brent2013,\n" \
"  author    = {Brent, R.P.},\n" \
"  title     = {Algorithms for Minimization Without Derivatives},\n" \
"  isbn      = {9780486143682},\n" \
"  series    = {Dover Books on Mathematics},\n" \
"  year      = {2013},\n" \
"  publisher = {Dover Publications}\n" \
"}\n"

class SchafferF6 : public nonlinearSystem {
public:

  SchafferF6()
  : nonlinearSystem( "Schaffer Function F6", SHAFFER_BIBTEX, 2 )
  { }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);

    real_type r = hypot(x1,x2);

    if ( r == 0.0 ) {
      f(0) = f(1) = 0;
      return;
    }

    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;

    real_type a  = 1/power2( 1.0 + 0.001 * r*r );
    real_type ar = - 0.004 * r / power3( 1.0 + 0.001 * r*r );

    real_type b  = power2(sin(r)) - 0.5;
    real_type br = sin ( 2.0 * r );

    f(0) = ( ar * b + a * br ) * rx1;
    f(1) = ( ar * b + a * br ) * rx2;

  }

  int_type
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {

    real_type x1 = x(0);
    real_type x2 = x(1);

    real_type r = hypot(x1,x2);

    if ( r == 0.0 ) {
      jac.setZero();
      return;
    }

    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;
    real_type r2  = r*r;
    real_type r3  = r2*r;

    real_type rx1x1 = x2*x2 / r3;
    real_type rx1x2 = - x1 * x2 / r3;
    real_type rx2x1 = - x1 * x2 / r3;
    real_type rx2x2 = x1*x1 / r3;

    real_type a   = 1/power2( 1.0 + 0.001 * r2 );
    real_type ar  = - 0.004 * r /power3( 1.0 + 0.001 * r2 );
    real_type arr = - 0.004 / power3( 1.0 + 0.001 * r2 ) 
                  + 0.000024 * r / power4( 1.0 + 0.001 * r2 );

    real_type b   = power2(sin(r))- 0.5;
    real_type br  = sin( 2.0 * r );
    real_type brr = 2.0 * cos( 2.0 * r );

    jac(0) = ( arr*b + 2*ar*br + a*brr ) * rx1*rx1 + ( ar*b + a*br ) * rx1x1;
    jac(1) = ( arr*b + 2*ar*br + a*brr ) * rx1*rx2 + ( ar*b + a*br ) * rx1x2;
    jac(2) = ( arr*b + 2*ar*br + a*brr ) * rx2*rx1 + ( ar*b + a*br ) * rx2x1;
    jac(3) = ( arr*b + 2*ar*br + a*brr ) * rx2*rx2 + ( ar*b + a*br ) * rx2x2;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  int_type
  numExactSolution() const override { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -5;
    x(1) = 10;
  }

  int_type
  numInitialPoint() const override { return 1; }

  void
  checkIfAdmissible( dvec_t const & X ) const override {
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SchafferF7 : public nonlinearSystem {
public:

  SchafferF7()
  : nonlinearSystem( "Schaffer Function F7", SHAFFER_BIBTEX, 2 )
  { }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);

    real_type r = hypot(x1,x2);

    if ( r == 0.0 ) {
      f(0) = f(1) = 0;
      return;
    }

    real_type a = sqrt ( r );
    real_type ar = 0.5 / sqrt ( r );

    real_type b = 1.0 + power2( sin ( 50.0 * pow(r,0.2) ) );
    real_type br = 10.0 * sin ( 100.0 * pow(r,0.2)  ) * pow(r,-0.8);

    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;

    f(0) = ( ar * b + a * br ) * rx1;
    f(1) = ( ar * b + a * br ) * rx2;
  }

  int_type
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {

    real_type x1 = x(0);
    real_type x2 = x(1);

    real_type r = hypot(x1,x2);

    if ( r == 0.0 ) {
      jac.setZero();
      return;
    }

    real_type rx1 = x1 / r;
    real_type rx2 = x2 / r;
    
    real_type r3 = r*r*r;

    real_type rx1x1 = x2*x2 / r3;
    real_type rx1x2 = - x1 * x2 / r3;
    real_type rx2x1 = - x1 * x2 / r3;
    real_type rx2x2 = x1*x1 / r3;
  
    //  F = A * B
    //  dFdX1 = ( Ar * B + A * Br ) * Rx1
    //  d2FdX1dX1 = ( Arr * B + Ar * Br ) * Rx1^2 + ( Ar * B + A * Br ) * Rx1x1
    //  etc
    real_type a   = sqrt ( r );
    real_type ar  = 0.5 / sqrt ( r );
    real_type arr = - 0.25 / sqrt ( r3 );

    real_type b   = 1.0 + power2( sin ( 50.0 * pow(r,0.2)  ) );
    real_type br  = 10.0 * sin ( 100.0 * pow(r,0.2) ) *  pow(r,-0.8);
    real_type brr = 200.0 * cos ( 100.0 * pow(r,0.2) ) *  pow(r,-1.6)
                  - 10.0 * sin ( 100.0 * pow(r,0.2) ) * 0.8 * pow(r,-1.8);

    jac(0) = ( arr*b + 2*ar*br + a*brr ) * rx1*rx1 + ( ar*b + a*br ) * rx1x1;
    jac(1) = ( arr*b + 2*ar*br + a*brr ) * rx1*rx2 + ( ar*b + a*br ) * rx1x2;
    jac(2) = ( arr*b + 2*ar*br + a*brr ) * rx2*rx1 + ( ar*b + a*br ) * rx2x1;
    jac(3) = ( arr*b + 2*ar*br + a*brr ) * rx2*rx2 + ( ar*b + a*br ) * rx2x2;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -5;
    x(1) = 10;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & X ) const override
  { }

};
