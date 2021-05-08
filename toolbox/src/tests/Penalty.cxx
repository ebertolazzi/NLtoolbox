/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define PENALTY_FUNCTION_BIBTEX \
"@techreport{Raydan:2004,\n" \
"  author = {William La Cruz and Jose Mario Martinez and Marcos Raydan},\n" \
"  title  = {Spectral residual method without gradient\n" \
"             information for solving large-scale nonlinear\n" \
"             systems of equations: Theory and experiments},\n" \
"  number = {Technical Report RT-04-08},\n" \
"  year   = {2004}\n" \
"}\n\n" \
"@article{LaCruz:2003,\n" \
"  author    = {William {La Cruz}  and  Marcos Raydan},\n" \
"  title     = {Nonmonotone Spectral Methods for Large-Scale Nonlinear Systems},\n" \
"  journal   = {Optimization Methods and Software},\n" \
"  year      = {2003},\n" \
"  volume    = {18},\n" \
"  number    = {5},\n" \
"  pages     = {583--599},\n" \
"  publisher = {Taylor & Francis},\n" \
"  doi       = {10.1080/10556780310001610493},\n" \
"}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class PenaltyIfunction : public nonlinearSystem {
public:

  PenaltyIfunction( int_type neq)
  : nonlinearSystem( "Penalty I", PENALTY_FUNCTION_BIBTEX, neq )
  { checkMinEquations(n,2); }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    if ( k == n-1 ) {
      real_type sum = 0;
      for ( int_type i = 0; i < n; ++i ) sum += x(i)*x(i);
      return (sum/n-1)/4;
    } else {
      return sqrt(1e-5)*(x(k)-1);
    }
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type sum = 0;
    for ( int_type i = 0; i < n; ++i ) sum += x(i)*x(i);
    for ( int_type i = 0; i < n-1; ++i ) f(i) = sqrt(1e-5)*(x(i)-1);
    f(n-1) = (sum/n-1)/4;
  }

  int_type
  jacobianNnz() const override {
    return 2*n-1;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n-1; ++i ) { SETIJ(i,i); }
    for ( int_type i = 0; i < n;   ++i ) { SETIJ(n-1,i); }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n-1; ++i ) jac(kk++) = sqrt(1e-5);
    for ( int_type i = 0; i < n;   ++i ) jac(kk++) = 0.5*x(i)/n;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 1.0/3.0;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class PenaltyN1 : public nonlinearSystem {
  real_type epsilon;
public:

  PenaltyN1( int_type neq )
  : nonlinearSystem(
      "Penalty Function #1",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      neq
    )
  , epsilon(0.00001)
  { checkMinEquations(n,2); }

  real_type
  sum( dvec_t const & x ) const {
    real_type t1 = 0;
    for ( int_type i = 0; i < n; ++i ) t1 += x(i)*x(i);
    return 4*t1 - 1;
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type ap = 2*epsilon;
    real_type t1 = sum(x);
    return (ap+t1)*x(k)-ap;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type ap = 2*epsilon;
    real_type t1 = sum(x);
    for ( int_type i = 0; i < n; ++i )
      f(i) = (ap+t1)*x(i)-ap;
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
    real_type ap = 2*epsilon;
    real_type t1 = sum(x);
    int_type  kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = 0; j < n; ++j ) {
        jac(kk) = 8*x(j)*x(i);
        if ( i == j ) jac(kk) += ap + t1;
        ++kk;
      }
    }
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = i+1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class PenaltyN2 : public nonlinearSystem {
  real_type epsilon;
public:

  PenaltyN2( int_type neq )
  : nonlinearSystem(
      "Penalty Function #2",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n",
      neq
    )
  , epsilon(0.00001)
  { checkMinEquations(n,2); }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type ap = epsilon;

    real_type t1 = -1.0;
    for ( int_type j = 0; j < n; ++j )
      t1 += ( n - j ) * (x(j)*x(j));

    real_type d2 = 1.0;
    real_type th = 4.0 * t1;
    real_type s2 = 0.0;
    for ( int_type j = 0; j < n; ++j ) {
      f(j) = ( n - j) * x(j) * th;
      real_type s1 = exp( x(j)/10.0);
      if ( j > 0 ) {
        real_type s3 = s1 + s2 - d2 * (exp( 0.1 ) + 1.0 );
        f(j)   += ap * s1 * ( s3 + s1 - 1.0 / exp ( 0.1 ) ) / 5.0;
        f(j-1) += ap * s2 * s3 / 5.0;
      }
      s2 = s1;
      d2 = d2 * exp ( 0.1 );
    }
    f(0) += 2.0 * ( x(0) - 0.2 );
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0; // fortran storage
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    real_type ap = 2*epsilon;

    real_type t1 = -1.0;
    for ( int_type j = 0; j < n; ++j )
      t1 += ( n - j ) * (x(j)*x(j));

    real_type d1 = exp( 0.1 );
    real_type d2 = 1.0;
    real_type s2 = 0.0;
    real_type th = 4.0 * t1;

    for ( int_type j = 0; j < n; ++j ) {

      jac[caddr(j,j)] = 8.0 * power2( (n - j) * x(j) ) + ( n - j ) * th;

      real_type s1 = exp ( x(j) / 10.0 );

      if ( j > 0 ) {
        real_type s3 = s1 + s2 - d2 * ( d1 + 1.0 );
        jac[caddr(j,j)]     += ap * s1 * ( s3 + s1 - 1.0 / d1 + 2.0 * s1 ) / 50.0;
        jac[caddr(j-1,j-1)] += ap * s2 * ( s2 + s3 ) / 50.0;
        for ( int_type k = 0; k < j; ++k ) {
          jac[caddr(j,k)] = 8.0 * (n - j) * (n - k) * x(j)*x(k);
        }
        jac[caddr(j,j-1)] += ap * s1 * s2 / 50.0;
      }
      s2 = s1;
      d2 = d1 * d2;
    }

    jac[caddr(0,0)] += 2;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = i+1; j < n; ++j )
        jac[caddr(i,j)] = jac[caddr(j,i)];
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type i = 0; i < n; ++i ) x(i) = 0.5;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
