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

class Chandrasekhar : public nonlinearSystem {
  dvec_t mu;
  real_type const w;
public:

  Chandrasekhar( real_type c, integer neq )
  : nonlinearSystem(
      "Chandrasekhar function",
      "@book{Kelley:1995,\n"
      "  author    = {Kelley, C.},\n"
      "  title     = {Iterative Methods for Linear and Nonlinear Equations},\n"
      "  publisher = {Society for Industrial and Applied Mathematics},\n"
      "  year      = {1995},\n"
      "  doi       = {10.1137/1.9781611970944},\n"
      "}\n\n"
      "@book{chandrasekhar1960,\n"
      "  author    = {Chandrasekhar, S.},\n"
      "  title     = {Radiative Transfer},\n"
      "  year      = {1960},\n"
      "  series    = {Dover Books on Intermediate and Advanced Mathematics},\n"
      "  publisher = {Dover Publications},\n"
      "  isbn      = {9780486605906}\n"
      "}\n",
      neq
    )
  , w(c/(2*neq))
  {
    mu.resize(neq);
    for ( integer i = 0; i < neq; ++i ) mu(i) = i + 0.5;
  }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    real_type tmp = 0;
    for ( integer j = 0; j < n; ++j )
      tmp += mu(j)*x(j)/(mu(i)+mu(j));
    return x(i)-1/(1-w*tmp);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n; ++i )
      f(i) = evalFk(x,i);
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i )
      for ( integer j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      real_type tmp = 0;
      for ( integer j = 0; j < n; ++j )
        tmp += mu(j)*x(j)/(mu(i)+mu(j));
      tmp = -w/power2(1-w*tmp);
      for ( integer j = 0; j < n; ++j ) {
        jac(kk) = tmp*mu(j)/(mu(i)+mu(j));
        if ( i == j ) jac(kk) += 1;
        ++kk;
      }
    }
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.fill(10);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
