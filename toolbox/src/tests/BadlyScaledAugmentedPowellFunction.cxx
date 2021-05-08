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

class BadlyScaledAugmentedPowellFunction : public nonlinearSystem {

  real_type
  phi( real_type t ) const {
    if ( t <= -1 )     return t/2-2;
    else if ( t >= 2 ) return t/2+2;
    else               return (-1924+t*(4551+t*(888-t*592)))/1998;
  }

  real_type
  phi_1( real_type t ) const {
    if ( t <= -1 )     return 0.5;
    else if ( t >= 2 ) return 0.5;
    else               return (4551+t*(2*888-t*3*592))/1998;
  }

public:

  BadlyScaledAugmentedPowellFunction( int_type neq )
  : nonlinearSystem(
      "Badly scaled augmented Powell’s function",
      "@article{Gasparo:2000,\n"
      "  Author    = {Maria Grazia Gasparo},\n"
      "  Title     = {A nonmonotone hybrid method for nonlinear systems},\n"
      "  Journal   = {Optimization Methods and Software},\n"
      "  Number    = {2},\n"
      "  Pages     = {79--94},\n"
      "  Publisher = {Taylor & Francis},\n"
      "  Volume    = {13},\n"
      "  Year      = {2000},\n"
      "  Doi       = {10.1080/10556780008805776},\n"
      "}\n",
      neq
    )
  { checkThree(n,3); }

  virtual
  real_type
  evalFk( dvec_t const & X, int_type k ) const override {
    int_type k1 = k % 3;
    int_type k2 = k - k1;
    //real_type const * x = X + k2;
    dvec_t const & x = X.segment(k2,3);
    switch ( k1 ) {
      case 0: return 10000 * (x(0) * x(1)) - 1.0;
      case 1: return exp(-x(1)) + exp(-x(0)) - 1.0001;
      case 2: return phi(x(2));
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & X, dvec_t & F ) const override {
    for ( int_type k = 0; k < n; k += 3 ) {
      dvec_t const & x = X.segment(k,3);
      F(k+0) = 10000 * (x(0) * x(1)) - 1.0;
      F(k+1) = exp(-x(1)) + exp(-x(0)) - 1.0001;
      F(k+2) = phi(x(2));
    }
  }

  virtual
  int_type
  jacobianNnz() const override
  { return (n/3)*5; }

  virtual
  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    int_type kk   = 0;
    int_type nblk = 0;
    for ( int_type k = 0; k < n; k += 3, nblk += 3 ) {
      i(kk) = nblk+0; j(kk) = nblk+0; ++kk;
      i(kk) = nblk+0; j(kk) = nblk+1; ++kk;
      i(kk) = nblk+1; j(kk) = nblk+0; ++kk;
      i(kk) = nblk+1; j(kk) = nblk+1; ++kk;
      i(kk) = nblk+2; j(kk) = nblk+2; ++kk;
    }
  }

  virtual
  void
  jacobian( dvec_t const & X, dvec_t & vals ) const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n; k += 3 ) {
      dvec_t const & x = X.segment(k,3);
      vals(kk) = 10000 * x(1); ++kk;
      vals(kk) = 10000 * x(0); ++kk;
      vals(kk) = -exp(-x(0));  ++kk;
      vals(kk) = -exp(-x(1));  ++kk;
      vals(kk) = phi_1(x(2));  ++kk;
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; k += 3 ) {
      x(k+0) = 0.109815932969981745568376164563E-4;
      x(k+1) = 9.10614673986652401094671049032;
      x(k+2) = 0.3998810580736440979319618294548679254646;
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 2; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch (idx) {
    case 0:
      for ( int_type k = 0; k < n; k += 3 ) {
        x(k+0) = 0;
        x(k+1) = 1;
        x(k+2) = -4;
      }
      break;
    case 1:
      for ( int_type k = 0; k < n; k += 3 ) {
        x(k+0) = 1e-3;
        x(k+1) = 18;
        x(k+2) = 1;
      }
      break;
    }
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 2; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
