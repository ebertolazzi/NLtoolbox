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

class HanbookFunction : public nonlinearSystem {
  mutable real_type sum1;
  mutable real_type sum2;
public:

  HanbookFunction( int_type neq)
  : nonlinearSystem(
      "Hanbook Function",
      "@techreport{Raydan:2004,\n"
      "  author = {William La Cruz and Jose Mario Martinez and Marcos Raydan},\n"
      "  title  = {Spectral residual method without gradient\n"
      "             information for solving large-scale nonlinear\n"
      "             systems of equations: Theory and experiments},\n"
      "  number = {Technical Report RT-04-08},\n"
      "  year   = {2004}\n"
      "}\n\n"
      "@article{LaCruz:2003,\n"
      "  author    = {William {La Cruz}  and  Marcos Raydan},\n"
      "  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems},\n"
      "  journal   = {Optimization Methods and Software},\n"
      "  year      = {2003},\n"
      "  volume    = {18},\n"
      "  number    = {5},\n"
      "  pages     = {583--599},\n"
      "  publisher = {Taylor & Francis},\n"
      "  doi       = {10.1080/10556780310001610493},\n"
      "}\n",
      neq
    )
  { checkMinEquations(n,2); }

  void
  sum( dvec_t const & x ) const {
    sum1 = sum2 = 0;
    for ( int_type i = 0; i < n; ++i ) {
      real_type xm = x(i) - 1;
      sum1 += xm;
      sum2 += xm*xm;
    }
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    sum(x);
    real_type S12 = sin( sum1 + sum2 );
    real_type S1  = 2 * sin(sum1);
    return 0.05*(x(k)-1) + (2+4*(x(k)-1)) * S12 + S1;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    sum(x);
    real_type S12 = sin( sum1 + sum2 );
    real_type S1  = 2 * sin(sum1);
    for ( int_type i = 0; i < n; ++i )
      f(i) = 0.05*(x(i)-1) + (2+4*(x(i)-1)) * S12 + S1;
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
    sum(x);
    int_type kk = 0;
    real_type S12 = sin( sum1 + sum2 );
    real_type C12 = cos( sum1 + sum2 );
    real_type C1  = 2 * cos(sum1);
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = 0; j < n; ++j ) {
        jac(kk) = (2+4*(x(i)-1)) * (2*x(j)-1) * C12 + C1;
        if ( i == j ) jac(kk) += 0.05 + 4 * S12;
        ++kk;
      }
    }
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 5;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
