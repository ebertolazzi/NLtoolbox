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

class BrownFunction : public nonlinearSystem {

public:

  BrownFunction()
  : nonlinearSystem(
      "Brown function",
      "@article{Qi:2006,\n"
      "  author  = {Qi, H. and Sun, D.},\n"
      "  title   = {A Quadratically Convergent Newton Method for\n"
      "             Computing the Nearest Correlation Matrix},\n"
      "  journal = {SIAM Journal on Matrix Analysis and Applications},\n"
      "  volume  = {28},\n"
      "  number  = {2},\n"
      "  pages   = {360--385},\n"
      "  year    = {2006},\n"
      "  doi     = {10.1137/050624509},\n"
      "}\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
     case 0: return power2(x(0))-(x(1)+1);
     case 1: return power2(x(0)-2)+power2(x(1)-0.5)-1;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(x(0))-(x(1)+1);
    f(1) = power2(x(0)-2)+power2(x(1)-0.5)-1;
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
    jac(0) =  2*x(0);
    jac(1) = -1;
    jac(2) =  2*x(0)-4;
    jac(3) =  2*x(1)-1;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1.06735;
    x(1) = 0.139228;
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.1;
    x(1) = 2;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
