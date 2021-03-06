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

class SSTnonlinearityTerm : public nonlinearSystem {
  integer   const idx;
  real_type const SCALE;
  real_type sst1;

public:

  SSTnonlinearityTerm( integer idx_in )
  : nonlinearSystem(
      fmt::format( "SST nonlinearity term, N.{}", idx_in ),
      "@article{Sincovec:1975,\n"
      "  author  = {Sincovec, Richard F. and Madsen, Niel K.},\n"
      "  title   = {Software for Nonlinear Partial Differential Equations},\n"
      "  journal = {ACM Trans. Math. Softw.},\n"
      "  volume  = {1},\n"
      "  number  = {3},\n"
      "  year    = {1975},\n"
      "  pages   = {232--260},\n"
      "  doi     = {10.1145/355644.355649}\n"
      "}\n",
      4
    )
  , idx(idx_in)
  , SCALE(1e7)
  {
    UTILS_ASSERT(
      idx == 0 || idx == 1,
      "SSTnonlinearityTerm, idx = {} must be 0 or 1", idx_in
    );
    switch ( idx ) {
      case 0: sst1 = 360;  break;
      case 1: sst1 = 3250; break;
    }
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 4E5-272.443800016*x(0)+0.1e-3*x(1)+0.7e-2*x(3)-3.67E-16*x(0)*x(1)-4.13E-12*x(0)*x(3);
    f(1) = 272.4438*x(0)-0.100016e-3*x(1)+3.67E-16*x(0)*x(1)-3.57E-15*x(1)*x(2);
    f(2) = -1.6E-8*x(2)+0.7e-2*x(3)+4.1283E-12*x(0)*x(3)-3.57E-15*x(1)*x(2)+800.0+sst1;
    f(3) = -0.7000016e-2*x(3)+3.57E-15*x(1)*x(2)-4.1283E-12*x(0)*x(3)+800.0;
    for ( integer i = 0; i < n; ++i ) f(i) /= SCALE;
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
    jac(kk++) = -272.443800016-3.67E-16*x(1)-4.13E-12*x(3);
    jac(kk++) = 0.0001-3.67E-16*x(0);
    jac(kk++) = 0.0;
    jac(kk++) = 0.007-4.13E-12*x(0);
    
    jac(kk++) = 272.4438+3.67E-16*x(1);
    jac(kk++) = -0.100016e-3+3.67E-16*x(0)-3.57E-15*x(2);
    jac(kk++) = -3.57E-15*x(1);
    jac(kk++) = 0.0;
    
    jac(kk++) = 4.1283E-12*x(3);
    jac(kk++) = -3.57E-15*x(2);
    jac(kk++) = -1.6E-8-3.57E-15*x(1);
    jac(kk++) = 0.007+4.1283E-12*x(0);
    
    jac(kk++) = -4.1283E-12*x(3);
    jac(kk++) = 3.57E-15*x(2);
    jac(kk++) = 3.57E-15*x(1);
    jac(kk++) = -0.7000016e-2-4.1283E-12*x(0);

    for ( integer i = 0; i < n*n; ++i ) jac(i) /= SCALE;
  }

  void
  getExactSolution( dvec_t & x, integer ns ) const override {
    switch ( ns ) {
      case 0:
        x(0) = 1.264224341494800168746220557617102243320e6;
        x(1) = 8.501430437421119197727722658565357159210e11;
        x(2) = 8.547006912362613482744674003963465220779e10;
        x(3) = 3.702993087637386517255325996036534779221e10;
      break;
      case 1:
        x(0) = 1.167055115589686729437819269080948198165e6;
        x(1) = 3.069755937686399353641443136774453022013e11;
        x(2) = 2.621168343680655428727664834155029399049e11;
        x(3) = 4.100816563193445712723351658449706009508e10;
      break;
    }
  }

  integer
  numExactSolution() const override
  { return 2; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1E9;
    x(1) = 1E9;
    x(2) = 1E13;
    x(3) = 1E7;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( integer i = 0; i < n; ++i )
      UTILS_ASSERT0( x(i) > 0, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(real_max);
    L.setZero();
  }
};
