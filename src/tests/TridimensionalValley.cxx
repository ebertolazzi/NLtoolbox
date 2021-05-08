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

class TridimensionalValley : public nonlinearSystem {

  real_type const c1;
  real_type const c2;

public:

  TridimensionalValley()
  : nonlinearSystem(
      "Tridimensional valley.",
      "no doc",
      3
    )
  , c1(1.003344481605351)
  , c2(-3.344481605351171E-3)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type bf = (c2*power2(x(0))+c1)*x(0);
    switch ( k ) {
      case 0: return bf*exp(-power2(x(0))/100)-1;
      case 1: return 10*(sin(x(0))-x(1));
      case 2: return 10*(cos(x(0))-x(2));
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type bf = (c2*power2(x(0))+c1)*x(0);
    f(0) = bf*exp(-power2(x(0))/100)-1;
    f(1) = 10*(sin(x(0))-x(1));
    f(2) = 10*(cos(x(0))-x(2));
  }

  int_type
  jacobianNnz() const override {
    return 5;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(1,0);
    SETIJ(1,1);
    SETIJ(2,0);
    SETIJ(2,2);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type t1 = x(0)*x(0);
    jac(0) = (c1-(t1/50.0)*(c1+c2*(t1-150.0)))*exp(-t1/100.0);
    jac(1) =  10*cos(x(0));
    jac(2) = -10;
    jac(3) = -10*sin(x(0));
    jac(4) = -10;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1.0103301175891008618821430258424435903873866121054053;
    x(1) = 0.847007375051043571769939585744456415641478463861070185;
    x(2) = 0.531581138312055623979884869864864195697816223034820704;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -4;
    x(1) = 1;
    x(2) = 2;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
