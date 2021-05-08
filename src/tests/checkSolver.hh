
/*

http://www.math.umn.edu/~olver/

*/

#ifndef CHECK_SOLVER_HH
#define CHECK_SOLVER_HH

#define USE_QR

#include "Alglin.hh"
#include "Alglin++.hh"
#include "NonlinearSystem.hh"
#include "testsNonlin.hh"
#include <fstream>

#ifndef F77NAME
  #define F77NAME(X) X##_
#endif

namespace NLsystem {

  using namespace std;
  using namespace ::CommonLoad;
  using namespace nonlinProblemLoad;

  typedef NonlinearSystem::NonlinearSystemBase SYSTEM_FUNCTIONS;
  
  static
  inline
  real_type
  maxabs( real_type a, real_type b )
  { return std::max( std::abs(a), std::abs(b) ); }

  static
  inline
  real_type
  maxabs3( real_type a, real_type b, real_type c )
  { return std::max( std::abs(a), maxabs(b,c) ); }
  
  static
  inline
  real_type
  abssum( real_type a, real_type b )
  { return std::abs(a)+std::abs(b); }

  static
  inline
  real_type
  abssum3( real_type a, real_type b, real_type c )
  { return std::abs(a)+std::abs(b)+std::abs(c); }

  class the_problem : public SYSTEM_FUNCTIONS {

    nonlinearSystem * p_problem;

    bool have_exact;
    bool use_qr;

    sizeType     n;
    real_type * x_stored;
    real_type * x_exact;
    real_type * JF;

    // soluzione minimi quadrati
    alglin::LU<real_type>            lu;
    alglin::QR<real_type>            qr;
    alglin::SVD<real_type>           svd;
    alglin::TridiagonalLU<real_type> tridLU;
    alglin::TridiagonalQR<real_type> tridQR;

    the_problem ( the_problem const & );
    the_problem const & operator = ( the_problem const & );

  public:

    the_problem()
    : SYSTEM_FUNCTIONS("solver")
    , p_problem(NULL)
    , x_stored(NULL)
    , x_exact(NULL)
    , JF(NULL)
    , use_qr(true)
    { }

    the_problem( nonlinearSystem * p )
    : SYSTEM_FUNCTIONS("solver")
    , p_problem(NULL)
    , x_stored(NULL)
    , x_exact(NULL)
    , JF(NULL)
    , use_qr(true)
    { setup( p ); }

    virtual
    ~the_problem() {
      delete [] x_stored;
      delete [] x_exact;
      delete [] JF;
    }

    void
    setup( nonlinearSystem * p ) {
      p_problem = p;
      n         = numEquations();
      try {
        if ( x_stored == NULL ) delete [] x_stored;
        if ( x_exact  == NULL ) delete [] x_exact;
        if ( JF       == NULL ) delete [] JF;

        x_stored = new real_type[n];
        x_exact  = new real_type[n];
        if ( p_problem->jac_tridiag() ) JF = new real_type[4*n];
        else                            JF = new real_type[n*n];
      } catch ( runtime_error & err ) {
        cerr << err.what() << '\n';
        exit(1);
      }
      have_exact = p_problem -> numExactSolution() > 0;
      if ( have_exact ) p_problem -> getExactSolution( x_exact, 0 );
    }

    bool
    jac_tridiag() const { return p_problem->jac_tridiag(); }

    virtual
    real_type
    evaluateComponent( dvec_t const & x, indexType k ) const {
      real_type tmp = p_problem -> evaluateComponent( x, k );
      checkNaN( &tmp, "evaluateComponent", 1, __LINE__, __FILE__ );
      return tmp;
    }

    virtual
    void
    evaluateFunction( dvec_t const & x, dvec_t & F ) const {
      p_problem -> evaluateFunction( x, F );
      checkNaN( F, "evaluateFunction", n, __LINE__, __FILE__ );
    }

    void
    evaluateJacobian( dvec_t const & x, real_type LU[] ) const {
      p_problem -> evaluateJacobian( x, LU );
      if ( p_problem->jac_tridiag() ) checkNaN( LU, "loadJacobian", 3*n, __LINE__, __FILE__ );
      else                            checkNaN( LU, "loadJacobian", n*n, __LINE__, __FILE__ );
    }

    virtual
    void
    loadJacobian( dvec_t const & x ) {
      p_problem -> evaluateJacobian( x, JF );
      if ( p_problem->jac_tridiag() ) {
        checkNaN( JF, "loadJacobian", 3*n, __LINE__, __FILE__ );
      } else {
        checkNaN( JF, "loadJacobian", n*n, __LINE__, __FILE__ );
      }
    }

    virtual
    void
    getJacobianScaling( real_type Df[], real_type Dx[] ) const {
      if ( p_problem->jac_tridiag() ) {
        real_type * L = JF+1;
        real_type * D = JF+n;
        real_type * U = JF+2*n;
        // calcolo bilanciatura
        Dx(0)   = 1/maxabs( D[0], L[0] );
        Dx[n-1] = 1/maxabs( D[n-1], U[n-2] );
        for ( indexType k = 1; k < n-1; ++k )
          Dx(k) = 1/maxabs3( L[k], D[k], U[k-1] );
        // calcolo bilanciatura
        Df(0)   = 1/maxabs( D[0]*Dx(0), U[0]*Dx(1) );
        Df(n-1) = 1/maxabs( D[n-1]*Dx[n-2], L[n-2]*Dx[n-1] );
        for ( indexType k = 1; k < n-1; ++k ) {
          Df(k) = 1/maxabs3( L[k-1]*Dx[k-1], D[k]*Dx(k), U[k]*Dx(k+1) );
        }
      } else {
        #if 1
        // calcolo bilanciatura
        //alglin::fill( n, Dx, 1, 1 );
        //alglin::fill( n, Df, 1, 1 );
        //alglin::equilibrate( n, n, JF, n, Df, Dx, 10, 1e-8 );
        for ( indexType k = 0; k < n; ++k )
          Dx(k) = sqrt(n/std::max(n*machineEps,alglin::asum( n, JF+k*n, 1 )));
        for ( indexType i = 0; i < n; ++i ) {
          Df(i) = 0;
          for ( indexType j = 0; j < n; ++j ) {
            real_type tmp = std::abs(JF[i+j*n]*Dx(j));
            if ( Df(i) < tmp ) Df(i) = tmp;
          }
          Df(i) = 1/std::max(sqrtMachineEps,Df(i));
        }
        //alglin::fill( n, Df, 1, 1 );
        #else
        // calcolo bilanciatura
        for ( indexType k = 0; k < n; ++k )
          Df(k) = (n/(sqrtMachineEps+alglin::asum( n, JF+k, n )));
        for ( indexType j = 0; j < n; ++j ) {
          Dx(j) = sqrtMachineEps;
          for ( indexType i = 0; i < n; ++i ) {
            real_type tmp = std::abs(JF[i+j*n]*Df(i));
            Dx(j) += tmp;
          }
          Dx(j) = sqrt(n/Dx(j));
        }
        #endif
      }
      //alglin::fill( n, Dx, 1, 1 );
      //alglin::fill( n, Df, 1, 1 );
    }

    virtual
    void
    scaleJacobian( real_type const Df[], real_type const Dx[] ) {
      if ( p_problem->jac_tridiag() ) {
        real_type * L = JF+1;
        real_type * D = JF+n;
        real_type * U = JF+2*n;
        // applico bilanciatura
        for ( indexType k = 0; k < n; ++k ) D[k] *= Dx(k)*Df(k);
        for ( indexType k = 0; k < n-1; ++k ) {
          L[k] *= Dx(k)*Df[k+1];
          U[k] *= Dx(k+1)*Df(k);
        }
      } else {
        // applico bilanciatura
        for ( indexType k = 0; k < n; ++k ) {
          alglin::scal( n, Df(k), JF+k,   n );
          alglin::scal( n, Dx(k), JF+k*n, 1 );
        }
      }
    }

    virtual
    void
    factorizeJacobian() {
      if ( p_problem->jac_tridiag() ) {
        real_type * L = JF+1;
        real_type * D = JF+n;
        real_type * U = JF+2*n;
        //tridLU.factorize( n, L, D, U );
        tridQR.factorize( n, L, D, U );
      } else {
        indexType info = lu.factorize( n, JF, n );
        NONLIN_ASSERT( info >= 0, "factorizeJacobian, getrf INFO = " << info );
        use_qr = info > 0 || lu.condInf( alglin::normInf( n, n, JF, n ) ) < 1e-6;
        //use_qr = false;
        cout << "COND = " << lu.condInf( alglin::normInf( n, n, JF, n ) ) << " info = " << info << '\n';
        if ( use_qr ) {
          //svd.factorize( n, n, JF_LU, n );
          qr.factorize( n, n, JF, n );
          qr.evalNumericRank( 1e-15 );
          cout << "NULLSPACE = " << n-qr.numericRank() << " n = " << n << '\n';
        }
      }
    }

/*\

  Sistema scalato  Df^(-1) F( Dx * y ) = 0
  
  ~~      -1                   ~~   -1
  JF  = Df   JF Dx --> JF = Df JF Dx

                             -1          ~~       -1
  direzione di ricerca d = JF   w = Dx ( JF )   Df  w

\*/

    virtual
    void
    solveJacobianSystem( real_type in_out[] ) const {
      if ( p_problem->jac_tridiag() ) {
        //tridLU.solve( in_out );
        tridQR.solve( in_out );
        //tridQR.minq( 1, in_out, n, 1e-15 );
      } else {
        if ( use_qr ) {
          //ls.solve( in_out, in_out );
          static int i = 1;
          i++;
          qr.minq( in_out, in_out, 0.0/(1+i*i*i*i) );
          //indexType rank = svd.solve( 1.0/(1+i*i*i*i), in_out, in_out );
          //cout << "NULLSPACE = " << n-rank << '\n';
        } else {
          lu.solve( in_out );
        }
      }
      checkNaN( in_out, "solveJacobianSystem", n, __LINE__, __FILE__ );
    }

/*\

  Sistema scalato  Df^(-1) F( Dx * y ) = 0
  
  ~~      -1                   ~~   -1
  JF  = Df   JF Dx --> JF = Df JF Dx

                                    -1   ~~
  direzione di ricerca d = JF w = Df   ( JF )  Dx  w

\*/

    virtual
    void
    multiplyByJacobian( real_type const in[], real_type out[] ) const {
      NONLIN_ASSERT( false, "NOT implemented" );
    }

    virtual
    sizeType
    evaluateMeritFunctions ( dvec_t const & x, real_type m[] ) const {
      NONLIN_ASSERT( false, "NOT implemented" );
    }

    virtual
    void
    initialPoint( dvec_t & x ) const {
      p_problem -> getInitialPoint( x, 0 ); // DA CAMBIARE @@@@@@@@@@@@@
    }

    virtual
    void
    checkIfAdmissible( dvec_t const & x ) const
    { p_problem -> checkIfAdmissible( x ); }

    virtual
    sizeType
    numEquations (void) const
    { return p_problem -> n_eqns(); }
  
    void
    check_solution(
      ostream      & s,
      dvec_t const & x,
      dvec_t const & f
    ) const {
      if ( have_exact ) {
        real_type maxnorm = 0;
        for ( indexType i = 0; i < numEquations(); ++i ) {
          real_type bf = x(i) - x_exact[i];
          if ( bf < 0       ) bf      = -bf;
          if ( bf > maxnorm ) maxnorm =  bf;
          if ( numEquations() <= 10 ) s << "\nxe[" << setw(2) << i << "]= " << x_exact[i];
        }
        s << "\n\n|| exact - numeric || = " << maxnorm << "\n\n";
      }

      if ( numEquations() <= 30 ) {
        for ( indexType i = 0; i < numEquations(); ++i )
          s << "x[" << setw(1) << i << "] = " << setprecision(18) << x(i) << "\n";
        for ( indexType i = 0; i < numEquations(); ++i )
          cout << "f[" << setw(1) << i << "] = " << setprecision(18) << f(i) << "\n";
      } else {
        real_type xinf = 0;
        real_type finf = 0;
        for ( indexType i = 0; i < numEquations(); ++i ) {
          xinf = std::max( xinf, std::abs(x(i)) );
          finf = std::max( finf, std::abs(f(i)) );
        }
        s << "||x||_inf = " << xinf << '\n';
        s << "||f||_inf = " << finf << '\n';
      }

    }

    string const & title() const { return p_problem -> title(); }
  };
}

#endif
