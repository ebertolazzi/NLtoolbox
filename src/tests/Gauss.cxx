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

class Gauss : public nonlinearSystem {
  real_type y[15];
public:

  Gauss()
  : nonlinearSystem(
      "Gaussian function",
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
    y[0]  = 0.0009;
    y[1]  = 0.0044;
    y[2]  = 0.0175;
    y[3]  = 0.0540;
    y[4]  = 0.1295;
    y[5]  = 0.2420;
    y[6]  = 0.3521;
    y[7]  = 0.3989;
    y[8]  = 0.3521;
    y[9]  = 0.2420;
    y[10] = 0.1295;
    y[11] = 0.0540;
    y[12] = 0.0175;
    y[13] = 0.0044;
    y[14] = 0.0009;
  }

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
    real_type x3 = x(2);
    f(0) = f(1) = f(2) = 0;
    for ( int_type i = 0; i < 15; ++i ) {
      real_type d1  = 0.5 * i;
      real_type d2  = 3.5 - d1 - x3;
      real_type arg = - 0.5 * x2 * d2 * d2;
      real_type t   = x1 * exp(arg) - y[i];

      f(0) += 2.0 * exp(arg) * t;
      f(1) -= x1 * exp(arg) * t * d2 * d2;
      f(2) += 2.0 * x1 * x2 * exp(arg) * t * d2;
    }

  }

  int_type
  jacobianNnz() const override {
    return 9;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I-1; jj(kk) = J-1; ++kk
    SETIJ(1,1);
    SETIJ(2,2);
    SETIJ(3,3);
    SETIJ(2,1);
    SETIJ(3,1);
    SETIJ(3,2);

    SETIJ(1,2);
    SETIJ(1,3);
    SETIJ(2,3);

    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type x1 = x(0);
    real_type x2 = x(1);
    real_type x3 = x(2);
    
    jac(0) = jac(1) = jac(2) = jac(3) = jac(4) = jac(5) = 0;

    for ( int_type i = 0; i < 15; ++i ) {
      real_type d1  = 0.5 * i;
      real_type d2  = 3.5 - d1 - x3;
      real_type arg = 0.5 * x2 * d2 * d2;
      real_type r   = exp(-arg);
      real_type t   = x1 * r - y[i];
      real_type t1  = 2.0 * x1 * r - y[i];
      
      real_type d2_2 = d2*d2;
      real_type d2_4 = d2_2*d2_2;

      jac(0) += r * r;
      jac(1) += r * t1 * d2_4;
      jac(2) += r * ( x2 * t1 * d2_2 - t );
      jac(3) -= r * t1 * d2_2;
      jac(4) += d2 * r * t1;
      jac(5) += d2 * r * ( t - arg * t1 );

    }

    jac(0) *= 2.0;
    jac(1) *= 0.5 * x1;
    jac(2) *= 2.0 * x1 * x2;
    jac(4) *= 2.0 * x2;
    jac(5) *= 2.0 * x1;

    jac(6) = jac(3);
    jac(7) = jac(4);
    jac(8) = jac(5);
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0.4;
    x(1) = 1;
    x(2) = 0;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 0;
    x(2) = 0;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( abs(x(i)) < 200, "x[" << i << "] = "<< x(i) << " too big");
  }

};
