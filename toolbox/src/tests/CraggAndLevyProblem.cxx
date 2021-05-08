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

class CraggAndLevyProblem : public nonlinearSystem {
public:
  
  CraggAndLevyProblem()
  : nonlinearSystem(
      "Cragg and Levy Problem",
      "@article{More:1981,\n"
      "  author  = {Mor{\'e}, Jorge J. and Garbow, Burton S. and Hillstrom, Kenneth E.},\n"
      "  title   = {Testing Unconstrained Optimization Software},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  year    = {1981},\n"
      "  volume  = {7},\n"
      "  number  = {1},\n"
      "  pages   = {17--41},\n"
      "  doi     = {10.1145/355934.355936},\n"
      "}\n",
      4
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return power2(exp(x(0))-x(1));
      case 1: return 10*power3(x(1)-x(2));
      case 2: return power2( tan( x(2)-x(3) ) );
      case 3: return x(3) - 1;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(exp(x(0))-x(1));
    f(1) = 10*power3(x(1)-x(2));
    f(2) = power2( tan( x(2)-x(3) ) );
    f(3) = x(3) - 1;
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {

    jac(0) = 2*(exp(x(0))-x(1))*exp(x(0));
    jac(1) = 2*(x(1)-exp(x(0)));
    jac(2) = 0;
    jac(3) = 0;

    jac(4) = 0;
    jac(5) = 30*power2(x(1)-x(2));
    jac(6) = -jac(0);
    jac(7) = 0;

    jac(8) = 0;
    jac(9) = 0;
    jac(10) = 2*sin(x(2)-x(3))/power3(cos(x(2)-x(3)));
    jac(11) = -jac(caddr(2,2));
    
    jac(12) = 0;
    jac(13) = 0;
    jac(14) = 0;
    jac(15) = 1;
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 4;
    x(1) = 2;
    x(2) = 2;
    x(3) = 2;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
