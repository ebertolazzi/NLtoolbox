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

class InfRefluxFunction : public nonlinearSystem {
public:

  InfRefluxFunction()
  : nonlinearSystem(
      "InfReflux function",
      "@article{Paterson:1986,\n"
      "  author  = {W.R. Paterson},\n"
      "  title   = {A new method for solving a class of nonlinear equations},\n"
      "  journal = {Chemical Engineering Science},\n"
      "  year    = {1986},\n"
      "  volume  = {41},\n"
      "  number  = {7},\n"
      "  pages   = {1935--1937},\n"
      "  doi     = {10.1016/0009-2509(86)87077-4}\n"
      "}\n",
      1
    )
  { }

  real_type
  evalFk( dvec_t const & x_in, int_type k ) const override {
    real_type x    = x_in[0];
    real_type arg1 = 1./(1.-x);
    real_type arg2 = 0.95-x;
    if ( x > 0 && arg1 > 0 && arg2 > 0 )
      return (1./63.)*log(x)+(64./63.)*log(arg1)+log(arg2)-log(0.9);
    return nan("InfRefluxFunction");
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x    = x_in[0];
    real_type arg1 = 1./(1.-x);
    real_type arg2 = 0.95-x;
    if ( x > 0 && arg1 > 0 && arg2 > 0 )
      f(0) = (1./63.)*log(x)+(64./63.)*log(arg1)+log(arg2)-log(0.9);
    else
      f(0) = nan("InfRefluxFunction");
  }

  int_type
  jacobianNnz() const override
  { return 1; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
  }

  void
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x    = x_in[0];
    real_type arg1 = 1./(1.-x);
    real_type arg2 = 0.95-x;
    if ( x > 0 && arg1 > 0 && arg2 > 0 )
      jac(0) = (1./63.)/x+(64./63.)/(1.-x)-1./(0.95-x);
    else
      jac(0) = nan("InfRefluxFunction");
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
      case 0 : x(0) = 0.23;  break; // default initial guess (critical as close to J=0 for x=0.229)
      case 1 : x(0) = 0.228; break; // an even harder initial point
      case 2 : x(0) = 0.6;   break; // and two relatively easy points
      case 3 : x(0) = 0.01;  break;
    }
  }

  int_type
  numInitialPoint() const override
  { return 4; }

  void
  checkIfAdmissible( dvec_t const & x_in ) const override {
    real_type x = x_in[0];
    NONLIN_ASSERT( x > 0 && x < 0.95, "ARGUMENT ERROR" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U[0] = 0.95; L[0] = 0;
  }

};
