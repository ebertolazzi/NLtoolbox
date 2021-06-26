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

class GeometricProgrammingFunction : public nonlinearSystem {
public:
  
  GeometricProgrammingFunction( integer  neq )
  : nonlinearSystem(
      "Geometric Programming Function",
      "@techreport{Raydan:2004,\n"
      "  author = {William La Cruz and Jose Mario Martinez and Marcos Raydan},\n"
      "  title  = {Spectral residual method without gradient\n"
      "             information for solving large-scale nonlinear\n"
      "             systems of equations: Theory and experiments},\n"
      "  number = {Technical Report RT-04-08},\n"
      "  year   = {2004}\n"
      "}\n",
      neq
    )
  {
    checkMinEquations(n,2);
  }

  real_type
  evalFk( dvec_t const & x, integer jj ) const override {
    real_type f = -1;
    for ( integer t = 1; t < 5; ++t ) {
      real_type t1 = 0.2*t;
      real_type t2 = t1-1;
      real_type tmp = 1;
      for ( integer k = 0; k < n; ++k ) {
        if ( jj != k ) tmp *= pow(x(k),t1);
        else           tmp *= t1*pow(x(k),t2);
      }
      f += tmp;
    }
    // t1 == 1
    // t2 == 0
    real_type tmp = 1;
    for ( integer k = 0; k < n; ++k ) {
      if ( jj != k ) tmp *= x(k);
    }
    f += tmp;
    return f;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f.fill(-1);
    for ( integer t = 1; t < 5; ++t ) {
      real_type t1 = 0.2*t;
      real_type t2 = t1-1;
      for ( integer i = 0; i < n; ++i ) {
        real_type tmp = 1;
        for ( integer k = 0; k < n; ++k ) {
          if ( i != k ) tmp *= pow(x(k),t1);
          else          tmp *= t1*pow(x(k),t2);
        }
        f(i) += tmp;
      }
    }
    for ( integer i = 0; i < n; ++i ) {
      real_type tmp = 1;
      for ( integer k = 0; k < n; ++k ) {
        if ( i != k ) tmp *= x(k);
      }
      f(i) += tmp;
    }
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
    jac.setZero();
    for ( integer t = 1; t < 5; ++t ) {
      real_type t1 = 0.2*t;
      real_type t2 = t1-1;
      real_type t3 = t2-1;
      integer kk = 0;
      for ( integer i = 0; i < n; ++i ) {
        for ( integer j = 0; j < n; ++j ) {
          real_type tmp = 1;
          if ( i == j ) {
            tmp *= t1*t2*pow(x(i),t3);
            for ( integer  k = 0; k < n; ++k )
              if ( i != k )
                tmp *= pow(x(k),t1);
          } else {
            for ( integer k = 0; k < n; ++k ) {
              if ( k == i || k == j ) {
                tmp *= t1*pow(x(k),t2);
              } else {
                tmp *= pow(x(k),t1);
              }
            }
          }
          jac(kk++) += tmp;
        }
      }
    }
    integer kk = 0;
    for ( integer i = 0; i < n; ++i ) {
      for ( integer j = 0; j < n; ++j ) {
        if ( i == j ) {
          ++kk;
        } else {
          real_type tmp = 1;
          for ( integer k = 0; k < n; ++k ) {
            if ( k == i || k == j ) continue;
            tmp *= x(k);
          }
          jac(kk++) += tmp;
        }
      }
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x.setZero();
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.fill(1);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < n-1; ++i )
      UTILS_ASSERT( x(i)>0, "x[{}] = {} must be > 0", i, x(i) );
  }

};
