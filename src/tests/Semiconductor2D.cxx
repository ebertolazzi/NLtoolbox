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

class Semiconductor2D : public nonlinearSystem {

  real_type alpha;
  real_type ni;
  real_type V;
  real_type D;

public:

  Semiconductor2D()
  : nonlinearSystem(
      "2D semiconductor",
      "@techreport{Nowak:1991,\n"
      "  author = {U. Nowak and L. Weimann},\n"
      "  title  = {A Family of Newton Co des for Systems of Highly Nonlinear Equations},\n"
      "  number = {Technical Report TR-91-10 (December 1991)},\n"
      "  year   = {1991}\n"
      "}\n",
      6
    )
  , alpha(38.683)
  , ni(1.22E10)
  , V(100)
  , D(1E7)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return exp(alpha*(x(2)-x(0))) - exp(alpha*(x(0)-x(1))) - D/ni;
      case 1: return x(1);
      case 2: return x(2);
      case 3: return exp(alpha*(x(5)-x(3))) - exp(alpha*(x(3)-x(4))) + D/ni;
      case 4: return x(4) - V;
      case 5: return x(5) - V;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = exp(alpha*(x(2)-x(0))) - exp(alpha*(x(0)-x(1))) - D/ni;
    f(1) = x(1);
    f(2) = x(2);
    f(3) = exp(alpha*(x(5)-x(3))) - exp(alpha*(x(3)-x(4))) + D/ni;
    f(4) = x(4) - V;
    f(5) = x(5) - V;
  }

  int_type
  jacobianNnz() const override {
    return 10;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);

    SETIJ(0,1);
    SETIJ(0,2);

    SETIJ(1,1);
    SETIJ(2,2);

    SETIJ(3,3);

    SETIJ(3,4);
    SETIJ(3,5);

    SETIJ(4,4);
    SETIJ(5,5);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = - alpha*( exp(alpha*(x(0)-x(1))) + exp(alpha*(x(2)-x(0))) );
    jac(1) = alpha*exp(alpha*(x(0)-x(1)));
    jac(2) = alpha*exp(alpha*(x(2)-x(0)));

    jac(3) = 1;
    jac(4) = 1;

    jac(5) = - alpha*( exp(alpha*(x(5)-x(3))) + exp(alpha*(x(3)-x(4))) );
    jac(6) = alpha*exp(alpha*(x(3)-x(4)));
    jac(7) = alpha*exp(alpha*(x(5)-x(3)));

    jac(8) = 1;
    jac(9) = 1;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
