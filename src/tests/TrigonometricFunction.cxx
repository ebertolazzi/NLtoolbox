/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

// TEST 220
/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class TrigonometricFunction : public nonlinearSystem {
public:
  
  TrigonometricFunction( int_type neq )
  : nonlinearSystem(
      "Trigonometric function",
      "Spedicato, E.\n"
      "Computational experience with quasi-newton algoritms.\n"
      "for minimization problems of moderately large size.\n"
      "Rep. CISE-N-175, Segrate (Milano), 1975.\n\n"
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
      neq
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
    real_type c_sum = 0;
    for ( int_type k = 0; k < n; ++k ) c_sum += cos(x(k));
    return n + (i+1) * (1-cos(x(i))) - sin(x(i)) - c_sum;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type c_sum = 0;
    for ( int_type i = 0; i < n; ++i ) c_sum += cos(x(i));
    for ( int_type i = 0; i < n; ++i ) {
      real_type t1 = n + (i+1) * (1-cos(x(i))) - sin(x(i)) - c_sum;
      real_type t2 = 2*sin(x(i)) - cos(x(i));
      f(i) = t1*t2;
    }
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
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      real_type c_sum = 0;
      for ( int_type j = 0; j < n; ++j ) c_sum += cos(x(j));
      real_type t1   = n + (i+1) * (1-cos(x(i))) - sin(x(i)) - c_sum;
      real_type t2   = 2*sin(x(i)) - cos(x(i));
      real_type t2_D = 2*cos(x(i)) + sin(x(i));
      for ( int_type j = 0; j < n; ++j ) {
        jac(kk) = t2*sin(x(j));
        if ( i == j )
          jac(kk) += t1*t2_D + t2*( (i+1)*sin(x(i)) - cos(x(i)) );
        ++kk;
      }
    }
  }

  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 1:
      x(0) = 0.9272952180016122324285124629224288040571;
      break;
    case 2:
      x(0) = 0.2430642022015621601022305882114319166998;
      x(1) = 0.6126761171373341819109378923948284268791;
      break;
    case 3:
      x(0) = 0.1386586620895952548146908578480828567861;
      x(1) = 0.1523812304815232863265755985362427737938;
      x(2) = 0.4677872324751889220566183484237867528208;
      break;
    case 4:
      x(0) = 0.8918060157524354158025129022089306168185e-1;
      x(1) = 0.9406975483676078600641830110116749612994e-1;
      x(2) = 0.1003490382176274133617689339270703727279;
      x(3) = 0.3808809535377544791201056909371821869283;
      break;
    case 5:
      x(0) = 0.6175491892349230228810295448928096965843e-1;
      x(1) = 0.6393992541842797329197797360800517769138e-1;
      x(2) = 0.6648669739226701899975920126598738991514e-1;
      x(3) = 0.6953055812954323530954539703632692011629e-1;
      x(4) = 0.3214890943456119400418258781551987027807;
      break;
    case 6:
      x(0) = 0.8132473891069603437697201795975588604775e-1;
      x(1) = 0.8530626471926446687284354166824938025959e-1;
      x(2) = 0.9026634644303587807017746608650930337132e-1;
      x(3) = 0.9681157292880582269931805627675287360652e-1;
      x(4) = 0.2883772090157352444000677796595946026999;
      x(5) = 0.2050122782185276818395942743633900930887;
      break;
    case 7:
      x(0) = 0.6359323463117480599802577468484466344791e-1;
      x(1) = -6.217264595163630379660558017324158265004;
      x(2) = 25.20139044537194916081602821610694937520;
      x(3) = -6.211249092913161143844503944299933134998;
      x(4) = 6.359242333673542087936017626710132641735;
      x(5) = 0.2487285997274638495927700030181851144571;
      x(6) = 0.1938449788520964162432781487951685607177;
      break;
    }
  }

  int_type
  numExactSolution() const override {
    if ( n < 8 ) return 1;
    else         return 0;
  }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k ) x(k) = 101.0 / real_type(100*n);
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for ( int_type i = 0; i < n; ++i )
    // NONLIN_ASSERT( abs(x(i)) < 2*m_pi, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    //for ( int_type i = 0; i < n; ++i )
    //  { U[i] = 2*m_pi; L[i] = -2*m_pi; }
    U.fill(real_max);
    L.fill(-real_max);
  }
};
