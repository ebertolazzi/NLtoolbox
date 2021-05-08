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

class GheriMancino : public nonlinearSystem {
  real_type alpha, beta, gamma;

public:

  // sum log(xi-2)^2+log(xi-10)^2 - prod( xi) ^(1/5)
  GheriMancino( int_type neq )
  : nonlinearSystem(
      "Gheri-Mancino function",
      "@Article{Gheri1971,\n"
      "  author  = {Gheri, G. and Mancino, O. G.},\n"
      "  title   = {A significant example to test methods for\n"
      "             solving systems of nonlinear equations},\n"
      "  journal = {CALCOLO},\n"
      "  year    = {1971},\n"
      "  volume  = {8},\n"
      "  number  = {1},\n"
      "  pages   = {107--113},\n"
      "  doi     = {10.1007/BF02575578}\n"
      "}\n",
      neq
    )
  , alpha(7)
  , beta(17)
  , gamma(4)
  {
  }

  real_type
  zfun( int_type i, int_type j, dvec_t const & x ) const {
    return sqrt( x(j)*x(j)+(i+1.0)/(j+1.0) );
  }

  real_type
  zfun_1( int_type i, int_type j, dvec_t const & x ) const {
    return x(j)/sqrt( x(j)*x(j)+(i+1.0)/(j+1.0) );
  }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type f = beta*n*x(i) + pow(i+1-0.5*n,gamma);
    for ( int_type j = 0; j < n; ++j ) {
      if ( i != j ) {
        real_type zij = zfun(i,j,x);
        real_type lij = log(zij);
        real_type ss  = sin(lij);
        real_type cc  = cos(lij);
        f += zij*( pow( ss, alpha ) + pow( cc, alpha ) );
      }
    }
    return f;
  }
  
  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type i = 0; i < n; ++i ) {
      f(i) = beta*n*x(i) + pow(i+1-0.5*n,gamma);
      for ( int_type j = 0; j < n; ++j ) {
        if ( i != j ) {
          real_type zij = zfun(i,j,x);
          real_type lij = log(zij);
          real_type ss  = sin(lij);
          real_type cc  = cos(lij);
          f(i) += zij*( pow( ss, alpha ) + pow( cc, alpha ) );
        }
      }
    }
  }

  virtual
  int_type
  jacobianNnz() const override
  { return n*n; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = 0; j < n; ++j ) {
        if ( i != j ) {
          real_type zij   = zfun(i,j,x);
          real_type zij_1 = zfun_1(i,j,x);
          real_type lij   = log(zij);
          real_type ss    = sin(lij);
          real_type cc    = cos(lij);
          jac(kk) = zij_1*( pow( ss, alpha ) + pow( cc, alpha ) )
                  + alpha*(pow( ss, alpha-1 )*cc - pow( cc, alpha-1 )*ss)*x(j)/zij;
        } else {
          jac(kk) = beta*n;
        }
        ++kk;
      }
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    real_type c = beta*n-(alpha+1)*(n-1);
    real_type k = beta*n+(alpha+1)*(n-1);
    for ( int_type i = 0; i < n; ++i ) {
      x(i) = (i+0.5*n) * gamma;
      for ( int_type j = 0; j < n; ++j ) {
        if ( i != j ) {
          real_type zij = sqrt( (i+1.0)/(j+1.0) );
          x(i) += zij*(pow( sin(log(zij)), alpha ) +
                       pow( cos(log(zij)), alpha ));
        
        }
      }
      x(i) *= -(c+k)/(2*k);
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
