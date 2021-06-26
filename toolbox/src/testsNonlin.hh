/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  This program is free software; you can redistribute it and/or modify    |
 |  it under the terms of the GNU General Public License as published by    |
 |  the Free Software Foundation; either version 2, or (at your option)     |
 |  any later version.                                                      |
 |                                                                          |
 |  This program is distributed in the hope that it will be useful,         |
 |  but WITHOUT ANY WARRANTY; without even the implied warranty of          |
 |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           |
 |  GNU General Public License for more details.                            |
 |                                                                          |
 |  You should have received a copy of the GNU General Public License       |
 |  along with this program; if not, write to the Free Software             |
 |  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               |
 |                                                                          |
 |  Copyright (C) 2003                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                | 
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Meccanica e Strutturale                  |
 |      Universita` degli Studi di Trento                                   |
 |      Via Mesiano 77, I-38050 Trento, Italy                               |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

/*
  http://www.sfu.ca/~ssurjano/optimization.html
  http://www-optima.amp.i.kyoto-u.ac.jp/member/student/hedar/Hedar_files/TestGO_files/Page364.htm
*/

#ifndef TESTS_NONLIN_HH
#define TESTS_NONLIN_HH

#include <iostream>
#include <sstream>

#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <cmath>
#include <map>

#define DEBUG 1
#define EIGEN_NO_AUTOMATIC_RESIZING 1

#include <Utils.hh>
#include <Eigen/Core>
#include <Eigen/Dense>

namespace NLproblem {

  typedef std::basic_ostream<char> ostream_type;

  using std::string;
  using std::map;
  using std::vector;
  using std::pair;
  using std::max;
  using std::min;

  using std::numeric_limits;

  typedef double  real_type;
  typedef int32_t integer;

  //! `m_e` the value of \f$ e \f$.
  static real_type const m_e = 2.718281828459045235360287471352662497757;

  //! `m_pi` the value of \f$ \pi \f$.
  static real_type const m_pi = 3.141592653589793238462643383279502884197;

  //! `m_2pi` the value of \f$ 2\pi \f$.
  static real_type const m_2pi = 6.283185307179586476925286766559005768394;

  //! `m_pi_2` the value of \f$ \pi/2 \f$.
  static real_type const m_pi_2 = 1.570796326794896619231321691639751442098;

  //! `m_pi_4` the value of \f$ \pi/4 \f$.
  static real_type const m_pi_4 = 0.7853981633974483096156608458198757210492;

  //! `m_1_pi` the value of \f$ 1/\pi \f$.
  static real_type const m_1_pi = 0.3183098861837906715377675267450287240689;

  //! `m_2_pi` the value of \f$ 2/\pi \f$.
  static real_type const m_2_pi = 0.6366197723675813430755350534900574481378;

  //! `m_sqrtpi` the value of \f$ \sqrt{\pi} \f$.
  static real_type const m_sqrtpi = 1.772453850905516027298167483341145182798;

  //! `m_2_sqrtpi` the value of \f$ 2/\sqrt{\pi} \f$.
  static real_type const m_2_sqrtpi = 1.128379167095512573896158903121545171688;

  //! `m_sqrt2` the value of \f$ \sqrt{2} \f$.
  static real_type const m_sqrt2 = 1.414213562373095048801688724209698078570;

  //! `m_1_sqrt2` the value of \f$ 1/\sqrt{2} \f$.
  static real_type const m_1_sqrt2 = 0.7071067811865475244008443621048490392850;

  static real_type real_max = numeric_limits<real_type>::max();

  typedef Eigen::Matrix<real_type,Eigen::Dynamic,Eigen::Dynamic> dmat_t;
  typedef Eigen::Matrix<real_type,Eigen::Dynamic,1>              dvec_t;
  typedef Eigen::Matrix<integer,Eigen::Dynamic,1>               ivec_t;

  class nonlinearBase {

    string const _title;
    string const _bibtex;

    nonlinearBase();
    nonlinearBase( nonlinearBase const & );
    nonlinearBase const & operator = ( nonlinearBase const & );

  protected:

    // U T I L I T Y  F U N C T I O N S -----------------------------------------

    real_type power2( real_type a ) const { return a*a; }
    real_type power3( real_type a ) const { return a*a*a; }
    real_type power4( real_type a ) const { real_type a2 = a*a; return a2*a2; }
    real_type power5( real_type a ) const { real_type a2 = a*a; return a2*a2*a; }
    real_type power6( real_type a ) const { real_type a2 = a*a; return a2*a2*a2; }

    void
    checkMinEquations( integer i, integer i_min ) const {
      UTILS_ASSERT(
        i >= i_min, "checkMinEquations:: i = {} < {}", i, i_min
      );
    }

    void
    checkEven( integer i, integer i_min ) const {
      UTILS_ASSERT(
        (i % 2) == 0 && i >= i_min, "checkEven:: odd index i = {}", i
      );
    }

    void
    checkOdd( integer i, integer i_min ) const {
      UTILS_ASSERT(
        (i % 2) != 0 && i >= i_min, "checkOdd:: odd index i = {}", i
      );
    }

    void
    checkThree( integer i, integer i_min ) const {
      UTILS_ASSERT(
        (i % 3) == 0 && i >= i_min, "checkThree:: index i = {}", i
      );
    }

    void
    checkFour( integer i, integer i_min ) const {
      UTILS_ASSERT(
        (i % 4) == 0 && i >= i_min, "checkFour:: index i = {}", i
      );
    }

  public:

    explicit nonlinearBase( string const & t, string const & b )
    : _title( t )
    , _bibtex( b )
    {}

    virtual ~nonlinearBase() {}

    string const & bibtex() const { return _bibtex; }
    string const & title()  const { return _title; }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class multivariateFunction : public nonlinearBase {

    multivariateFunction( multivariateFunction const & );
    multivariateFunction const & operator = ( multivariateFunction const & );

  protected:
    integer n;

  public:

    multivariateFunction( string const & t, string const & b, integer _n )
    : nonlinearBase( fmt::format("{} neq = {}", t, _n ), b )
    , n(_n)
    { }

    virtual ~multivariateFunction() {}

    virtual real_type eval( dvec_t const & x ) const = 0;
    virtual void      gradient( dvec_t const & x, dvec_t & g ) const = 0;

    virtual integer hessianNnz() const = 0;
    virtual void    hessian( dvec_t const & x, dvec_t & jac ) const = 0;
    virtual void    hessianPattern( ivec_t & i, ivec_t & j ) const = 0;

    virtual integer numExactSolution() const = 0;
    virtual void    getExactSolution( dvec_t & x, integer idx ) const = 0;

    virtual integer numInitialPoint() const = 0;
    virtual void    getInitialPoint( dvec_t & x, integer idx ) const = 0;

    virtual void    checkIfAdmissible( dvec_t const & x ) const = 0;

    virtual void
    boundingBox( dvec_t & L, dvec_t & U ) const {
      L.fill( -real_max );
      U.fill( real_max );
    }

    integer dimX( void ) const { return n; }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class nonlinearLeastSquares: public nonlinearBase {

    nonlinearLeastSquares( nonlinearLeastSquares const & );
    nonlinearLeastSquares const & operator = ( nonlinearLeastSquares const & );

  protected:

    integer n, m;

  public:

    nonlinearLeastSquares(
      string const & t,
      string const & b,
      integer       _n,
      integer       _m
    )
    : nonlinearBase(
        fmt::format( "{} dimF = {}, dimX = {}", t, _n, _m ), b
      )
    , n(_n)
    , m(_m)
    {}

    virtual ~nonlinearLeastSquares() {}

    virtual real_type evalFk( dvec_t const & x, integer k ) const = 0;
    virtual void      evalF ( dvec_t const & x, dvec_t & f ) const = 0;

    virtual integer jacobianNnz() const = 0;
    virtual void    jacobian( dvec_t const & x, dvec_t & jac ) const = 0;
    virtual void    jacobianPattern( ivec_t & i, ivec_t & j ) const = 0;

    virtual integer tensorNnz() const = 0;
    virtual void    tensor( dvec_t const & x, dvec_t const & lambda, dvec_t & jac ) const = 0;
    virtual void    tensorPattern( ivec_t & i, ivec_t & j ) const = 0;

    virtual integer numExactSolution() const = 0;
    virtual void    getExactSolution( dvec_t & x, integer idx ) const = 0;

    virtual integer numInitialPoint() const = 0;
    virtual void    getInitialPoint( dvec_t & x, integer idx ) const = 0;

    virtual void checkIfAdmissible( dvec_t const & x ) const = 0;

    virtual void
    boundingBox( dvec_t & L, dvec_t & U ) const {
      L.fill( -real_max );
      U.fill( real_max );
    }

    integer dimF( void ) const { return n; }
    integer dimX( void ) const { return m; }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class nonlinearSystem: public nonlinearBase {
  
    nonlinearSystem( nonlinearSystem const & );
    nonlinearSystem const & operator = ( nonlinearSystem const & );

  protected:

    // U T I L I T Y  F U N C T I O N S -----------------------------------------
  
    // n = dimensione matrice
    // i = riga
    // j = colonna
    // indirizzamento fortran
    integer faddr( integer i, integer j ) const
    { return (i-1) + (j-1) * n; }
    
    integer caddr( integer i, integer j ) const
    { return i + j * n; }

    integer n;

  public:

    nonlinearSystem( string const & t, string const & b, integer _n )
    : nonlinearBase( fmt::format("{} neq = {}", t, _n ), b )
    , n(_n)
    { }

    virtual ~nonlinearSystem() {}

    //void
    //setup( string const & t, integer _n )
    //{ theTitle = t; n = _n; }

    virtual real_type evalFk( dvec_t const & x, integer k ) const = 0;
    virtual void      evalF ( dvec_t const & x, dvec_t & f ) const = 0;

    virtual integer jacobianNnz() const = 0;
    virtual void    jacobian( dvec_t const & x, dvec_t & jac ) const = 0;
    virtual void    jacobianPattern( ivec_t & i, ivec_t & j ) const = 0;

    virtual integer numExactSolution() const = 0;
    virtual void    getExactSolution( dvec_t & x, integer idx ) const = 0;

    virtual integer numInitialPoint() const = 0;
    virtual void    getInitialPoint( dvec_t & x, integer idx ) const = 0;

    virtual void checkIfAdmissible( dvec_t const & x ) const = 0;

    virtual void
    boundingBox( dvec_t & L, dvec_t & U ) const {
      L.fill( -real_max );
      U.fill( real_max );
    }

    integer numEqns( void ) const { return n; }

    integer
    fill_CSR(
      dvec_t const & x,
      ivec_t       & R,
      ivec_t       & J,
      dvec_t       & values
    ) const;

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  class nonlinearSystemFromMultivariateFunction: public nonlinearSystem {

    nonlinearSystemFromMultivariateFunction(
      nonlinearSystemFromMultivariateFunction const &
    );
    nonlinearSystemFromMultivariateFunction const &
    operator = (nonlinearSystemFromMultivariateFunction const &);

    multivariateFunction const * pMF;

  public:

    nonlinearSystemFromMultivariateFunction(
      multivariateFunction const * _pMF
    )
    : nonlinearSystem( _pMF->title(), _pMF->bibtex(), _pMF->dimX() )
    , pMF(_pMF) {
    }

    virtual ~nonlinearSystemFromMultivariateFunction() {}

    virtual
    real_type
    evalFk( dvec_t const & x, integer k ) const {
      dvec_t g(n);
      evalF( x, g );
      return g(k);
    }

    virtual
    void
    evalF( dvec_t const & x, dvec_t & g ) const {
      pMF->gradient( x, g );
    }

    virtual
    integer
    jacobianNnz() const
    { return pMF->hessianNnz(); }

    virtual
    void
    jacobian( dvec_t const & x, dvec_t & hess ) const
    { pMF->hessian(x,hess); }

    virtual
    void
    jacobianPattern( ivec_t & i, ivec_t & j ) const
    { pMF->hessianPattern( i, j ); }

    virtual
    integer
    numExactSolution() const
    { return pMF->numExactSolution(); }

    virtual
    void
    getExactSolution( dvec_t & x, integer idx ) const
    { pMF->getExactSolution( x, idx ); }

    virtual
    integer
    numInitialPoint() const
    { return pMF->numInitialPoint(); }

    virtual
    void
    getInitialPoint( dvec_t & x, integer idx ) const
    { pMF->getInitialPoint( x, idx ); }

    virtual
    void
    checkIfAdmissible( dvec_t const & x ) const
    { pMF->checkIfAdmissible( x ); }

    virtual
    void
    boundingBox( dvec_t & L, dvec_t & U ) const
    { pMF->boundingBox( L, U ); }

    integer numEqns( void ) const { return pMF->dimX(); }

    string const & title(void) const { return pMF->title(); }

  };

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#if 0
  class nonlinearSystemFromLeastSquares: public nonlinearSystem {

    nonlinearSystemFromLeastSquares(nonlinearSystemFromLeastSquares const &);
    nonlinearSystemFromLeastSquares const &
    operator = (nonlinearSystemFromLeastSquares const &);

    nonlinearLeastSquares const * pLS;

  public:

    nonlinearSystemFromLeastSquares( nonlinearLeastSquares const * _pLS )
    : nonlinearSystem( _pLS->title(), _pLS->dimX() )
    , pLS(_pLS) {
    }

    virtual ~nonlinearSystemFromLeastSquares() {}

    virtual
    real_type
    evalFk( dvec_t const & x, integer k ) const {
      dvec_t g(n);
      evalF( x, g );
      return g(k);
    }

    virtual
    void
    evalF( dvec_t const & x, dvec_t & g ) const {
      pLS->gradient( x, g );
    }

    virtual
    integer
    jacobianNnz() const
    { return pMF->hessianNnz(); }

    virtual
    void
    jacobian( dvec_t const & x, dmat_t & hess ) const
    { pMF->hessian(x,hess); }

    virtual
    void
    jacobianPattern( ivec_t & i, ivec_t & j ) const
    { pMF->hessianPattern(i,j); }

    virtual
    integer
    numExactSolution() const
    { return pMF->numExactSolution(); }

    virtual
    void
    getExactSolution( dvec_t & x, integer idx ) const
    { pMF->getExactSolution(x,idx); }

    virtual
    integer
    numInitialPoint() const
    { return pMF->numInitialPoint(); }

    virtual
    void
    getInitialPoint( dvec_t & x, integer idx ) const
    { pMF->getInitialPoint(x,idx); }

    virtual
    void
    checkIfAdmissible( dvec_t const & x ) const
    { pMF->checkIfAdmissible(x); }

    virtual
    void
    boundingBox( dvec_t & L, dvec_t & U ) const
    { pMF->boundingBox(L,U); }

    integer numEqns( void ) const { return pMF->dimX(); }

    string const & title(void) const { return pMF->title(); }

  };
#endif

  extern vector<nonlinearSystem*> theProblems;
  extern map<string,integer>      theProblemsMap;
  void initProblems();

}

#endif
