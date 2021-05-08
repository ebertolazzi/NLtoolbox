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

class BrownAndGearhartFunction : public nonlinearSystem {
  real_type const sqrt2;
public:

  BrownAndGearhartFunction()
  : nonlinearSystem(
      "Brown and Gearhart function",
      "@Article{Brow1971,\n"
      "  author  = {Brow, Kenneth M. and Gearhart, William B.},\n"
      "  title   = {Deflation techniques for the calculation of further solutions\n"
      "             of a nonlinear system},\n"
      "  journal = {Numerische Mathematik},\n"
      "  year    = {1971},\n"
      "  volume  = {16},\n"
      "  number  = {4},\n"
      "  pages   = {334--342},\n"
      "  doi     = \"10.1007/BF02165004\",\n"
      "}\n",
      3
    )
  , sqrt2(sqrt(2.0)) {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
     case 0: return power2(x(0))+2*power2(x(1))-4;
     case 1: return power2(x(0))+power2(x(1))+x(2)-8;
     case 2: return power2(x(0)-1)+power2(2*x(1)-sqrt2)+power2(x(2)-5)-4;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(x(0))+2*power2(x(1))-4;
    f(1) = power2(x(0))+power2(x(1))+x(2)-8;
    f(2) = power2(x(0)-1)+power2(2*x(1)-sqrt2)+power2(x(2)-5)-4;
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
    jac(0) = 2*x(0);
    jac(1) = 4*x(1);
    jac(2) = 0;

    jac(3) = 2*x(0);
    jac(4) = 2*x(1);
    jac(5) = 1;

    jac(6) = 2*(x(0)-1);
    jac(7) = 4*(2*x(1)-sqrt2);
    jac(8) = 2*(x(2)-5);
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = sqrt2;
    x(2) = 6;
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 0.7;
    x(2) = 5;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
