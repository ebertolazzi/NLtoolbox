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

class Shubert : public nonlinearSystem {
public:

  Shubert()
  : nonlinearSystem(
      "Shubert Function",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      2
    )
  {
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = f(1) = 0;

    real_type factor1 = 0.0;
    real_type df1dx1  = 0.0;
    for ( int_type i = 0; i < 5; ++i ) {
      real_type y = i+1;
      factor1 += y * cos ( ( y + 1.0 ) * x(0)+ y );
      df1dx1  -= y * ( y + 1.0 ) * sin ( ( y + 1.0 ) * x(0)+ y );
    }

    real_type factor2 = 0.0;
    real_type df2dx2  = 0.0;
    for ( int_type i = 0; i < 5; ++i ) {
      real_type y = i;
      factor2 += y * cos ( ( y + 1.0 ) * x(1) + y );
      df2dx2  -= y * ( y + 1.0 ) * sin ( ( y + 1.0 ) * x(1) + y );
    }
    f(0) = df1dx1 * factor2;
    f(1) = factor1 * df2dx2;
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
    real_type factor1   = 0.0;
    real_type df1dx1    = 0.0;
    real_type factor1_D = 0.0;
    real_type df1dx1_D  = 0.0;
    for ( int_type i = 0; i < 5; ++i ) {
      real_type y  = i+1;
      real_type y1 = y + 1.0;
      factor1   += y * cos ( y1 * x(0)+ y );
      df1dx1    -= y * y1 * sin ( y1 * x(0)+ y );

      factor1_D -= y * sin( y1 * x(0)+ y ) * y1;
      df1dx1_D  -= y * power2(y1) * cos ( y1 * x(0)+ y );
    }

    real_type factor2   = 0.0;
    real_type df2dx2    = 0.0;
    real_type factor2_D = 0.0;
    real_type df2dx2_D  = 0.0;
    for ( int_type i = 0; i < 5; ++i ) {
      real_type y  = i;
      real_type y1 = y + 1.0;
      factor2   += y * cos ( y1 * x(1) + y );
      df2dx2    -= y * y1 * sin ( y1 * x(1) + y );
      factor2_D -= y * sin ( y1 * x(1) + y ) * y1;
      df2dx2_D  -= y * power2(y1) * cos ( y1 * x(1) + y );
    }
    jac(0) = df1dx1_D * factor2;
    jac(1) = df1dx1 * factor2_D;
    jac(2) = factor1_D * df2dx2;
    jac(3) = factor1 * df2dx2_D;
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
    x(0) = 0.5;
    x(1) = 1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
