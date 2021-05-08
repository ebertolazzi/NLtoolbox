/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

/*
"@Article{Spedicato1997,\n"
"  author  = {Spedicato, E. and Huang, Z.},\n"
"  title   = {Numerical experience with newton-like methods\n"
"             for nonlinear algebraic systems},\n"
"  journal = {Computing},\n"
"  year    = {1997},\n"
"  volume  = {58},\n"
"  number  = {1},\n"
"  pages   = {69--89},\n"
"  doi     = {10.1007/BF02684472},\n"
"}\n"

  Numerical Experience with Newton-like Methods for Nonlinear Algebraic Systems
  E. Spedicato, Z. Huang
  COmputing N. 58, 1997
  
  Test Problems for Unconstrained Optimization
  Ladislav Luksan Jan Vlcek
  Technical report No. 897
  January 2003

  Sparse Test Problems for Unconstrained Optimization
  Ladislav Luksan Jan Vlcek
  Technical report No. 1064
  January 2010
 
  A MODIFIED NEWTON METHOD FOR SOLVING NON-LINEAR ALGEBRAIC EQUATIONS
  Satya N. Atluri*, Chein-Shan Liu**, and Chung-Lun Kuo***
  Journal of Marine Science and Technology, Vol. 17, No. 3, pp. 238-247 (2009)
 */

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

#define RKM_BIBTEX \
"@book{meresoo:1990,\n" \
"  title     = {Test Examples of Systems of Nonlinear Equations: Version 3-90},\n" \
"  author    = {Meresoo, T. and Roose, A. and Kulla,\n" \
"               V. and Estonian Software and Computer Service Company},\n" \
"  year      = 1990,\n" \
"  publisher = {Estonian Software and Computer Service Company}\n" \
"}\n"

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo129 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo129()
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.129",RKM_BIBTEX,3)
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
    case 0: return x(0) + x(1) - 2;
    case 1: return x(0) - log(x(1)) + x(2) - 2;
    case 2: return x(1)*x(1) - 2*x(2) + 1;
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0) + x(1) - 2;
    f(1) = x(0) - log(x(1)) + x(2) - 2;
    f(2) = x(1)*x(1) - 2*x(2) + 1;
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
    jac(0) = 1;
    jac(1) = 1;
    jac(2) = 0;

    jac(3) = 1;
    jac(4) = -1/x(1);
    jac(5) = 1;

    jac(6) = 0;
    jac(7) = 2*x(1);
    jac(8) = -2;
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
      case 0:
        x(0) = x(1) = x(2) = 1;
      break;
      case 1:
        x(0) = -0.285888702509893511;
        x(1) = 2.28588870250989329;
        x(2) = 3.11264358013118247;
      break;
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 2; }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(0.5);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

  string note() const { return "Jacobian is singular at the solution"; }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo201 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo201( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.201",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    if ( k == 0 ) f = 1-x(0);
    else          f = 10*k*power2(x(k)-x(k-1));
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 1-x(0);
    for ( int_type k = 1; k < n; ++k )
    	f(k) = 10*k*power2(x(k)-x(k-1));
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 2*n-1;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    for ( int_type k = 1; k < n; ++k ) {
      SETIJ(k,k-1);
      SETIJ(k,k);
    }

    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = -1;
    for ( int_type k = 1; k < n; ++k ) {
      jac[2*k-1] = -20*k*(x(k)-x(k-1));
      jac[2*k]   =  20*k*(x(k)-x(k-1));
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.fill(1);
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(-1.2);
    x(n-1) = -1;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

  string note() const { return "Jacobian is singular at the solution"; }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo202 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo202( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.202",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    if ( k == n-1 ) f = x(n-1)-0.1*power2(x(0));
    else            f = x(k)-0.1*power2(x(k+1));
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n-1; ++k )
    	f(k) = x(k) - 0.1*power2(x(k+1));
    f(n-1) = x(n-1) - 0.1*power2(x(0));
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 2*n;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 0; k < n-1; ++k ) {
      SETIJ(k,k);
      SETIJ(k,k+1);
    }
    SETIJ(n-1,n-1);
    SETIJ(n-1,0);
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n-1; ++k ) {
      jac(kk++) = 1;
      jac(kk++) = -0.2*x(k+1);
    }
    jac(kk++) = 1;
    jac(kk++) = -0.2*x(0);
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
      case 0: x.setZero(); break;
      case 1: x.fill(10);  break;
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 2; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(2);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo203 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo203( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.203",RKM_BIBTEX,neq)
  { checkMinEquations(n,2); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    for ( int_type i = 0; i < n; ++i ) f += x(i);
    if ( k == n-1 ) f -= 1;
    else            f += x(k)-(n+1);
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
  	real_type sum = 0;
    for ( int_type i = 0; i < n;   ++i ) sum += x(i);
    for ( int_type k = 0; k < n-1; ++k ) f(k) = x(k)-(n+1)+sum;
    f(n-1) = sum - 1;
  }

  virtual
  int_type
  jacobianNnz() const override
  { return n*n; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0; // fortran storing
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.fill(1);
    for ( int_type i = 0; i < n-1; ++i ) jac[caddr(i,i)] += 1;
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2: return 2;
    case 5: return 3;
    case 10:
    case 20:
    case 30: return 1;
    default: return 0;
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      switch ( idx ) {
        case 0: x(0) = x(1) = 1;      break;
        case 1: x(0) = 0.5; x(1) = 2; break;
      }
      break;
    case 5:
      switch ( idx ) {
        case 0:
          x(0) = x(1) = x(2) = x(3) = x(4) = 1;
          break;
        case 1:
         x(0) = x(1) = x(2) = x(3) = 0.91635458253384926;
         x(4) = 1.4182270873307543;
         break;
        case 2:
         x(0) = x(1) = x(2) = x(3) = -5.7904308849411582E-1;
         x(4) = 9.8952154424705789;
         break;
      }
      break;
    case 10:
    case 20:
    case 30:
      x.fill(1);
    }

  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(0.5);
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo204 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo204( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.204",RKM_BIBTEX,neq)
  { checkEven(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    if ( (k % 2) == 0 ) f = 1-x(k);
    else                f = 10*(x(k)-power2(x(k-1)));
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; k += 2 ) f(k) = 1-x(k);
    for ( int_type k = 1; k < n; k += 2 ) f(k) = 10*(x(k)-power2(x(k-1)));
  }

  virtual
  int_type
  jacobianNnz() const override {
    return n+n/2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 0; k < n; k += 2 ) {
      SETIJ(k,k);
    }
    for ( int_type k = 1; k < n; k += 2 ) {
      SETIJ(k,k-1);
      SETIJ(k,k);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n; k += 2 ) jac(kk++) = -1;
    for ( int_type k = 1; k < n; k += 2 ) {
      jac(kk++) = -20*x(k-1);
      jac(kk++) = 10;
    }
  }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.fill(1);
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; k += 2 ) x(k) = -1.2;
    for ( int_type k = 1; k < n; k += 2 ) x(k) = 1;
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo205 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo205( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.205",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    for ( int_type j = 0; j < n; ++j ) f += power3(x(j));
    f = x(k) - 0.5*(f+k+1)/n;
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
  	real_type acc = 0;
    for ( int_type j = 0; j < n; ++j ) acc += power3(x(j));
    acc /= 2*n;
    for ( int_type j = 0; j < n; ++j ) f(j) = x(j) - acc - (0.5/n)*(j+1);
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
    real_type bf = -1.5/n;
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = 0; j < n; ++j ) {
        jac(kk) = bf*power2(x(j));
        if ( i == j ) jac(kk) += 1;
        ++kk;
      }
    }
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 3;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      switch ( idx ) {
      case 0:
        x(0) = 1.01242808774686965;
        x(1) = 1.26242808774686965;
        break;
      case 1:
        x(0) = -1.68508582524678752;
        x(1) = -1.43508582524678752;
        break;
      case 2:
        x(0) = 0.297657737499917863;
        x(1) = 0.547657737499917863;
        break;
      }
      break;
    case 5:
      switch ( idx ) {
      case 0:
        x(0) = 1;
        x(1) = 1.1;
        x(2) = 1.2;
        x(3) = 1.3;
        x(4) = 1.4;
        break;
      case 1:
        x(0) = -1.72736184954957039;
        x(1) = -1.62736184954957053;
        x(2) = -1.52736184954957044;
        x(3) = -1.42736184954957057;
        x(4) = -1.32736184954957048;
        break;
      case 2:
        x(0) = 0.127361849549570388;
        x(1) = 0.227361849549570394;
        x(2) = 0.327361849549570427;
        x(3) = 0.427361849549570405;
        x(4) = 0.527361849549570438;
        break;
      }
      break;
    case 10:
      switch ( idx ) {
      case 0:
        x(0) = 0.994471118343573157;
        x(1) = 1.04447111834357309;
        x(2) = 1.09447111834357313;
        x(3) = 1.14447111834357318;
        x(4) = 1.194471118343573;
        x(5) = 1.24447111834357305;
        x(6) = 1.29447111834357309;
        x(7) = 1.34447111834357313;
        x(8) = 1.39447111834357318;
        x(9) = 1.444471118343573;
        break;
      case 1:
        x(0) = -1.74181474184792218;
        x(1) = -1.69181474184792235;
        x(2) = -1.64181474184792231;
        x(3) = -1.59181474184792227;
        x(4) = -1.54181474184792222;
        x(5) = -1.49181474184792218;
        x(6) = -1.44181474184792213;
        x(7) = -1.39181474184792231;
        x(8) = -1.34181474184792227;
        x(9) = -1.29181474184792222;
        break;
      case 2:
        x(0) = 0.0723436235043492942;
        x(1) = 0.122343623504349297;
        x(2) = 0.172343623504349314;
        x(3) = 0.222343623504349303;
        x(4) = 0.272343623504349264;
        x(5) = 0.322343623504349364;
        x(6) = 0.372343623504349297;
        x(7) = 0.422343623504349341;
        x(8) = 0.472343623504349275;
        x(9) = 0.522343623504349264;
        break;
      }
      break;
    case 20:
      switch ( idx ) {
      case 0:
        x(0) = 0.991518329652893882;
        x(1) = 1.0165183296528939;
        x(2) = 1.04151832965289382;
        x(3) = 1.06651832965289395;
        x(4) = 1.09151832965289386;
        x(5) = 1.11651832965289377;
        x(6) = 1.1415183296528939;
        x(7) = 1.16651832965289382;
        x(8) = 1.19151832965289395;
        x(9) = 1.21651832965289386;
        x(10) = 1.24151832965289377;
        x(11) = 1.2665183296528939;
        x(12) = 1.29151832965289382;
        x(13) = 1.31651832965289395;
        x(14) = 1.34151832965289386;
        x(15) = 1.36651832965289377;
        x(16) = 1.3915183296528939;
        x(17) = 1.41651832965289382;
        x(18) = 1.44151832965289395;
        x(19) = 1.46651832965289386;
        break;
      case 1:
        x(0) = -1.74911100354538962;
        x(1) = -1.72411100354538971;
        x(2) = -1.69911100354538958;
        x(3) = -1.67411100354538966;
        x(4) = -1.64911100354538975;
        x(5) = -1.62411100354538962;
        x(6) = -1.59911100354538971;
        x(7) = -1.57411100354538958;
        x(8) = -1.54911100354538966;
        x(9) = -1.52411100354538975;
        x(10) = -1.49911100354538962;
        x(11) = -1.47411100354538971;
        x(12) = -1.44911100354538958;
        x(13) = -1.42411100354538966;
        x(14) = -1.39911100354538975;
        x(15) = -1.37411100354538962;
        x(16) = -1.34911100354538971;
        x(17) = -1.32411100354538958;
        x(18) = -1.29911100354538966;
        x(19) = -1.27411100354538975;
        break;
      case 2:
        x(0) = 0.0450926738924959589;
        x(1) = 0.0700926738924959603;
        x(2) = 0.0950926738924959686;
        x(3) = 0.120092673892495963;
        x(4) = 0.145092673892495971;
        x(5) = 0.170092673892495994;
        x(6) = 0.19509267389249596;
        x(7) = 0.220092673892495982;
        x(8) = 0.245092673892495949;
        x(9) = 0.270092673892495971;
        x(10) = 0.295092673892495994;
        x(11) = 0.320092673892496016;
        x(12) = 0.345092673892495982;
        x(13) = 0.370092673892496005;
        x(14) = 0.395092673892495971;
        x(15) = 0.420092673892495994;
        x(16) = 0.445092673892496016;
        x(17) = 0.470092673892495982;
        x(18) = 0.495092673892496005;
        x(19) = 0.520092673892496027;
        break;
      }
      break;
    case 30:
      switch ( idx ) {
      case 0:
        x(0) = 0.990509030546815938;
        x(1) = 1.00717569721348266;
        x(2) = 1.02384236388014926;
        x(3) = 1.04050903054681587;
        x(4) = 1.0571756972134827;
        x(5) = 1.07384236388014931;
        x(6) = 1.09050903054681592;
        x(7) = 1.10717569721348252;
        x(8) = 1.12384236388014935;
        x(9) = 1.14050903054681596;
        x(10) = 1.15717569721348257;
        x(11) = 1.1738423638801494;
        x(12) = 1.19050903054681601;
        x(13) = 1.20717569721348261;
        x(14) = 1.22384236388014922;
        x(15) = 1.24050903054681605;
        x(16) = 1.25717569721348266;
        x(17) = 1.27384236388014926;
        x(18) = 1.29050903054681587;
        x(19) = 1.3071756972134827;
        x(20) = 1.32384236388014931;
        x(21) = 1.34050903054681592;
        x(22) = 1.35717569721348252;
        x(23) = 1.37384236388014935;
        x(24) = 1.39050903054681596;
        x(25) = 1.40717569721348257;
        x(26) = 1.4238423638801494;
        x(27) = 1.44050903054681601;
        x(28) = 1.45717569721348261;
        x(29) = 1.47384236388014922;
        break;
      case 1:
        x(0) = -1.75155356027278342;
        x(1) = -1.73488689360611681;
        x(2) = -1.71822022693944998;
        x(3) = -1.70155356027278337;
        x(4) = -1.68488689360611676;
        x(5) = -1.66822022693945016;
        x(6) = -1.65155356027278333;
        x(7) = -1.63488689360611672;
        x(8) = -1.61822022693945011;
        x(9) = -1.60155356027278351;
        x(10) = -1.58488689360611668;
        x(11) = -1.56822022693945007;
        x(12) = -1.55155356027278346;
        x(13) = -1.53488689360611685;
        x(14) = -1.51822022693945002;
        x(15) = -1.50155356027278342;
        x(16) = -1.48488689360611681;
        x(17) = -1.4682202269394502;
        x(18) = -1.45155356027278337;
        x(19) = -1.43488689360611676;
        x(20) = -1.41822022693945016;
        x(21) = -1.40155356027278355;
        x(22) = -1.38488689360611672;
        x(23) = -1.36822022693945011;
        x(24) = -1.35155356027278351;
        x(25) = -1.33488689360611668;
        x(26) = -1.31822022693945007;
        x(27) = -1.30155356027278346;
        x(28) = -1.28488689360611685;
        x(29) = -1.26822022693945002;
        break;
      case 2:
        x(0) = 0.0360445297259673891;
        x(1) = 0.0527111963926340521;
        x(2) = 0.0693778630593007289;
        x(3) = 0.0860445297259673919;
        x(4) = 0.102711196392634055;
        x(5) = 0.119377863059300732;
        x(6) = 0.136044529725967395;
        x(7) = 0.152711196392634058;
        x(8) = 0.169377863059300721;
        x(9) = 0.186044529725967384;
        x(10) = 0.202711196392634047;
        x(11) = 0.219377863059300737;
        x(12) = 0.2360445297259674;
        x(13) = 0.252711196392634063;
        x(14) = 0.269377863059300726;
        x(15) = 0.286044529725967389;
        x(16) = 0.302711196392634052;
        x(17) = 0.319377863059300715;
        x(18) = 0.336044529725967378;
        x(19) = 0.352711196392634041;
        x(20) = 0.369377863059300704;
        x(21) = 0.386044529725967367;
        x(22) = 0.40271119639263403;
        x(23) = 0.419377863059300748;
        x(24) = 0.436044529725967411;
        x(25) = 0.452711196392634074;
        x(26) = 0.469377863059300737;
        x(27) = 0.4860445297259674;
        x(28) = 0.502711196392634063;
        x(29) = 0.519377863059300782;
        break;
      }
      break;
    }
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(1.5);
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo206 : public nonlinearSystem {
  real_type h2;
public:
  
  RooseKullaLombMeressoo206( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.206",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); h2 = 1.0/power2(n+1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type xc = x(k);
  	real_type xp = k < n-1 ? x(k+1) : 0;
  	real_type xm = k > 0   ? x(k-1) : 0;
    return xp-2*xc+xm - h2 *exp(xc);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k ) {
    	real_type xc = x(k);
    	real_type xp = k < n-1 ? x(k+1) : 0;
    	real_type xm = k > 0   ? x(k-1) : 0;
      f(k) = xp-2*xc+xm - h2 *exp(xc);
    }
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 1; k < n;   ++k ) { SETIJ(k,k-1); }
    for ( int_type k = 0; k < n-1; ++k ) { SETIJ(k,k+1); }
    for ( int_type k = 0; k < n;   ++k ) { SETIJ(k,k); }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 1; k < n;   ++k ) jac(kk++) = 1;
    for ( int_type k = 0; k < n-1; ++k ) jac(kk++) = 1;
    for ( int_type k = 0; k < n;   ++k ) jac(kk++) = -2 - h2 *exp(x(k));
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = -0.100488400337317069;
      x(1) = -0.100488400337317069;
      break;
    case 5:
      x(0) = -0.0635730237796020142;
      x(1) = -0.101079225590438845;
      x(2) = -0.113478165704209086;
      x(3) = -0.101079225590438845;
      x(4) = -0.0635730237796020142;
      break;
    case 10:
      x(0) = -0.0380470822832765579;
      x(1) = -0.068138233862133038;
      x(2) = -0.0905092917544646769;
      x(3) = -0.10533104513306768;
      x(4) = -0.112714562777264757;
      x(5) = -0.112714562777264757;
      x(6) = -0.10533104513306768;
      x(7) = -0.0905092917544646769;
      x(8) = -0.068138233862133038;
      x(9) = -0.038047082283276551;
      break;
    case 20:
      x(0) = -0.0209484001798467684;
      x(1) = -0.0396762346150324255;
      x(2) = -0.0562247027012899511;
      x(3) = -0.0706295728406963474;
      x(4) = -0.0829215019487556798;
      x(5) = -0.0931263031751204673;
      x(6) = -0.101265166939055559;
      x(7) = -0.107354839384380743;
      x(8) = -0.11140776149579644;
      x(9) = -0.113432171358331196;
      x(10) = -0.113432171358331182;
      x(11) = -0.111407761495796426;
      x(12) = -0.107354839384380729;
      x(13) = -0.101265166939055545;
      x(14) = -0.0931263031751204534;
      x(15) = -0.0829215019487556521;
      x(16) = -0.0706295728406963197;
      x(17) = -0.0562247027012899234;
      x(18) = -0.0396762346150324047;
      x(19) = -0.020948400179846758;
      break;
    case 30:
      x(0) = -0.0144369806182624797;
      x(1) = -0.0278482934603717126;
      x(2) = -0.0402476022487183935;
      x(3) = -0.0516473776603994031;
      x(4) = -0.0620589494503979225;
      x(5) = -0.0714925529960212941;
      x(6) = -0.0799573706823796804;
      x(7) = -0.087461568494601391;
      x(8) = -0.094012328134088427;
      x(9) = -0.0996158749325536247;
      x(10) = -0.104277501798100969;
      x(11) = -0.108001589391558919;
      x(12) = -0.110791622698076339;
      x(13) = -0.112650204128129788;
      x(14) = -0.113579063253107349;
      x(15) = -0.113579063253107335;
      x(16) = -0.112650204128129788;
      x(17) = -0.110791622698076339;
      x(18) = -0.108001589391558905;
      x(19) = -0.104277501798100955;
      x(20) = -0.0996158749325536108;
      x(21) = -0.0940123281340884132;
      x(22) = -0.0874615684946013633;
      x(23) = -0.0799573706823796526;
      x(24) = -0.0714925529960212663;
      x(25) = -0.0620589494503978947;
      x(26) = -0.0516473776603993753;
      x(27) = -0.0402476022487183727;
      x(28) = -0.0278482934603716953;
      x(29) = -0.014436980618262471;
      break;
    }
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.setZero();
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo207 : public nonlinearSystem {
  real_type gamma;
public:
  
  RooseKullaLombMeressoo207( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.207",RKM_BIBTEX,neq)
  , gamma(0.1)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type xc = x(k);
  	real_type xp = k < n-1 ? x(k+1) : 0;
  	real_type xm = k > 0   ? x(k-1) : 0;
    return (3-gamma*xc)*xc+1-xm-2*xp;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k ) {
    	real_type xc = x(k);
    	real_type xp = k < n-1 ? x(k+1) : 0;
    	real_type xm = k > 0   ? x(k-1) : 0;
      f(k) = (3-gamma*xc)*xc+1-xm-2*xp;
    }
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 1; k < n;   ++k ) { SETIJ(k,k-1); }
    for ( int_type k = 0; k < n-1; ++k ) { SETIJ(k,k+1); }
    for ( int_type k = 0; k < n;   ++k ) { SETIJ(k,k); }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 1; k < n;   ++k ) jac(kk++) = -1;
    for ( int_type k = 0; k < n-1; ++k ) jac(kk++) = -2;
    for ( int_type k = 0; k < n;   ++k ) jac(kk++) = 3-2*gamma*x(k);
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 2;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      switch ( idx ) {
      case 0:
        x(0) = 20.9686582909684525;
        x(1) = 9.96875591028267927;
        break;
      case 1:
        x(0) = -0.685453889794504945;
        x(1) = -0.551673186443478292;
        break;
      }
      break;
    case 5:
      switch ( idx ) {
      case 0:
        x(0) = 10.8597777768098975;
        x(1) = 10.8929279971301387;
        x(2) = 5.47670908975717818;
        x(3) = 1.76888251337925362;
        x(4) = 0.25852195788334581;
        break;
      case 1:
        x(0) = -1.5293511879989905;
        x(1) = -1.91097253481018181;
        x(2) = -1.78437400965571991;
        x(3) = -1.38027427739523056;
        x(4) = -0.773482265306932204;
        break;
      }
      break;
    case 10:
      switch ( idx ) {
      case 0:
        x(0) = 14.4359473997111909;
        x(1) = 11.7340922332053701;
        x(2) = 3.99871862308393045;
        x(3) = -0.16845571330670267;
        x(4) = -1.75346174786930264;
        x(5) = -2.19969617021264563;
        x(6) = -2.16474654344672679;
        x(7) = -1.88157810993199437;
        x(8) = -1.41701070236339088;
        x(9) = -0.785122965109708248;
        break;
      case 1:
        x(0) = -1.96909750815374274;
        x(1) = -2.64751351206147811;
        x(2) = -2.83718790384275144;
        x(3) = -2.8345068598189691;
        x(4) = -2.73488779472511823;
        x(5) = -2.55905882466501389;
        x(6) = -2.29858344303975581;
        x(7) = -1.93252004445795178;
        x(8) = -1.43622003127863795;
        x(9) = -0.791206423601281572;
        break;
      }
      break;
    case 20:
      switch ( idx ) {
      case 0:
        x(0) = 15.2764638454963269;
        x(1) = 11.7461783871066672;
        x(2) = 3.58240032282524856;
        x(3) = -0.641168312964383369;
        x(4) = -2.27350747113667895;
        x(5) = -2.84811886128854175;
        x(6) = -3.04101360876585014;
        x(7) = -3.09984917093945933;
        x(8) = -3.1097201961549672;
        x(9) = -3.0981736936814257;
        x(10) = -3.07233445425563545;
        x(11) = -3.03137678448305348;
        x(12) = -2.97036021007190376;
        x(13) = -2.88100391174524884;
        x(14) = -2.75133493955649211;
        x(15) = -2.56499265094333007;
        x(16) = -2.3007808716064142;
        x(17) = -1.93335461289545485;
        x(18) = -1.43653448650018212;
        x(19) = -0.79130598984776257;
        break;
      case 1:
        x(0) = -2.04511475622648753;
        x(1) = -2.77679685264649745;
        x(2) = -3.02816793889987679;
        x(3) = -3.11234353533562258;
        x(4) = -3.13876544765076693;
        x(5) = -3.14456883057665504;
        x(6) = -3.14188617855131014;
        x(7) = -3.13411729048722565;
        x(8) = -3.1213674059817329;
        x(9) = -3.10213918788524357;
        x(10) = -3.07368845588766515;
        x(11) = -3.03184112608173084;
        x(12) = -2.97052049186879019;
        x(13) = -2.88105977439294048;
        x(14) = -2.75135468683627105;
        x(15) = -2.56499977369673227;
        x(16) = -2.3007835090801767;
        x(17) = -1.93335561455466376;
        x(18) = -1.43653486390840945;
        x(19) = -0.79130610934649992;
        break;
      }
      break;
    case 30:
      switch ( idx ) {
      case 0:
        x(0) = 15.2931719204216012;
        x(1) = 11.745702511253814;
        x(2) = 3.57389143252621277;
        x(3) = -0.650649105411801876;
        x(4) = -2.28408658729946801;
        x(5) = -2.86165790515736784;
        x(6) = -3.05989786239379979;
        x(7) = -3.12716658742612319;
        x(8) = -3.14975949321800242;
        x(9) = -3.15710518936978879;
        x(10) = -3.15914369628296443;
        x(11) = -3.15917239442777253;
        x(12) = -3.15820525438590227;
        x(13) = -3.15643470580650387;
        x(14) = -3.1537034341177943;
        x(15) = -3.14963006579175797;
        x(16) = -3.14360185919570911;
        x(17) = -3.13469938835462081;
        x(18) = -3.12156516570161857;
        x(19) = -3.10220650856120583;
        x(20) = -3.07371144107897365;
        x(21) = -3.0318490084888472;
        x(22) = -2.97052321270752451;
        x(23) = -2.88106072267857449;
        x(24) = -2.75135502205215898;
        x(25) = -2.56499989460753319;
        x(26) = -2.30078355385205269;
        x(27) = -1.93335563155811641;
        x(28) = -1.43653487031502247;
        x(29) = -0.791306111375025267;
        break;
      case 1:
        x(0) = -2.04665682215098732;
        x(1) = -2.77942544060933994;
        x(2) = -3.03207003883383752;
        x(3) = -3.11806477396577719;
        x(4) = -3.14717853826406069;
        x(5) = -3.15697205799869751;
        x(6) = -3.1601924466152429;
        x(7) = -3.16114345590571721;
        x(8) = -3.16126035799173088;
        x(9) = -3.16099716158523725;
        x(10) = -3.16046071615948687;
        x(11) = -3.15961809036597874;
        x(12) = -3.15835610131762179;
        x(13) = -3.15648576992995533;
        x(14) = -3.15372072502463707;
        x(15) = -3.14963592314447371;
        x(16) = -3.14360384462249964;
        x(17) = -3.13470006195778073;
        x(18) = -3.12156539454732718;
        x(19) = -3.10220658646387104;
        x(20) = -3.07371146767713377;
        x(21) = -3.03184901761026149;
        x(22) = -2.97052321585604062;
        x(23) = -2.88106072377591715;
        x(24) = -2.75135502244006558;
        x(25) = -2.56499989474744927;
        x(26) = -2.30078355390386191;
        x(27) = -1.93335563157779244;
        x(28) = -1.4365348703224361;
        x(29) = -0.791306111377372723;
        break;
      }
      break;
    }
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(-1.0);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( int_type i = 3; i < n; ++i )
      NONLIN_ASSERT( x(i) < 0, "Bad range" );
  }

  virtual
  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    int_type i = 0;
    for (; i < 3; ++i )
      { U[i] = real_max; L[i] = -real_max; }
    for (; i < n; ++i )
      { U[i] = 0; L[i] = -real_max; }
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo208 : public nonlinearSystem {
  real_type const K1;
  real_type const K2;
  real_type const K3;
  int_type  const R1;
  int_type  const R2;
public:
  
  RooseKullaLombMeressoo208( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.208",RKM_BIBTEX,neq)
  , K1(1)
  , K2(1)
  , K3(1)
  , R1(3)
  , R2(3)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type sum = 0;
    for ( int_type i = std::max(0,k-R1); i <= std::min(n-1,k+R2); ++i )
      if ( i != k ) sum += x(i)*(1+x(i));
  	return (K1+K2*x(k)*x(k))*x(k)+1-K3*sum;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k )
      f(k) = evalFk( x, k );
  }

  virtual
  int_type
  jacobianNnz() const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      kk += std::min(n-1,i+R2) - std::max(0,i-R1) + 1;
    }
    return kk;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = std::max(0,i-R1); j <= std::min(n-1,i+R2); ++j ) {
        SETIJ(i,j);
      }
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type j = std::max(0,i-R1); j <= std::min(n-1,i+R2); ++j ) {
        if ( i == j ) jac(kk) = K1 + 3*K2*x(i)*x(i);
        else          jac(kk) = -K3*(1+2*x(j));
        ++kk;
      }
    }
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 1;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = -0.754877666246692725;
      x(1) = -0.754877666246692725;
      break;
    case 5:
      x(0) = -0.808025022386814507;
      x(1) = -0.87166824608729554;
      x(2) = -0.87166824608729554;
      x(3) = -0.87166824608729554;
      x(4) = -0.808025022386814507;
      break;
    case 10:
      x(0) = -0.801937946359539633;
      x(1) = -0.842558051527937502;
      x(2) = -0.881216883188142397;
      x(3) = -0.911897857643132359;
      x(4) = -0.891348792903065235;
      x(5) = -0.891348792903065235;
      x(6) = -0.911897857643132359;
      x(7) = -0.881216883188142397;
      x(8) = -0.842558051527937502;
      x(9) = -0.801937946359539633;
      break;
    case 20:
      x(0) = -0.800413820399077491;
      x(1) = -0.839125061489862767;
      x(2) = -0.880560005368091314;
      x(3) = -0.920665554037870759;
      x(4) = -0.897476024922029492;
      x(5) = -0.88270771158442396;
      x(6) = -0.880029836830879542;
      x(7) = -0.889821686957802283;
      x(8) = -0.891604540812876101;
      x(9) = -0.889489266395734002;
      x(10) = -0.889489266395734002;
      x(11) = -0.891604540812876101;
      x(12) = -0.889821686957802283;
      x(13) = -0.880029836830879542;
      x(14) = -0.88270771158442396;
      x(15) = -0.897476024922029492;
      x(16) = -0.920665554037870759;
      x(17) = -0.880560005368091314;
      x(18) = -0.839125061489862656;
      x(19) = -0.800413820399077491;
      break;
    case 30:
      x(0) = -0.800390085020429187;
      x(1) = -0.839050608964517775;
      x(2) = -0.880547772685447039;
      x(3) = -0.920819039192321442;
      x(4) = -0.897634652952583934;
      x(5) = -0.882511895121084633;
      x(6) = -0.879442500798295668;
      x(7) = -0.889801872319090914;
      x(8) = -0.892712487130870924;
      x(9) = -0.890760710534785227;
      x(10) = -0.887785410621690585;
      x(11) = -0.888000317543117035;
      x(12) = -0.889093931815407745;
      x(13) = -0.889479357496824941;
      x(14) = -0.889087453456860799;
      x(15) = -0.889087453456860799;
      x(16) = -0.88947935749682483;
      x(17) = -0.889093931815407634;
      x(18) = -0.888000317543117035;
      x(19) = -0.887785410621690585;
      x(20) = -0.890760710534785227;
      x(21) = -0.892712487130871035;
      x(22) = -0.889801872319090914;
      x(23) = -0.879442500798295668;
      x(24) = -0.882511895121084633;
      x(25) = -0.897634652952583822;
      x(26) = -0.920819039192321553;
      x(27) = -0.880547772685447039;
      x(28) = -0.839050608964517775;
      x(29) = -0.800390085020429187;
      break;
    }
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(-1);
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo209 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo209( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.209",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }
  
  void
  checkx( dvec_t const & x ) const {
    for ( int_type i = 0; i < n; ++i )
      NONLIN_ASSERT( x(i) > 0,
                     "RooseKullaLombMeressoo209, found x[" << i <<
                     "] = " << x(i) << " <= 0" );
  }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    if ( k == 0 ) f = power2(x(0))-1;
    else          f = power2(x(k-1))-1+log(x(k));
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(x(0))-1;
    for ( int_type k = 1; k < n; ++k )
      f(k) = power2(x(k-1))-1+log(x(k));
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 2*n-1;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    for ( int_type k = 1; k < n; ++k ) {
      SETIJ(k,k-1);
      SETIJ(k,k);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = 2*x(0);
    for ( int_type k = 1; k < n; ++k ) {
      jac(kk++) = 2*x(k-1);
      jac(kk++) = 1/x(k);
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.fill(1);
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(0.5);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( int_type i = 1; i < n; ++i ) {
      NONLIN_ASSERT( x(i) > 0, "Bad range" );
      NONLIN_ASSERT( std::abs(x(i)) < 1000, "Bad range" );
    }
  }

  virtual
  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(1000);
    L.fill(-1000);
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

// DA RICONTROLLARE!
class RooseKullaLombMeressoo210 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo210( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.210",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    for ( int_type j = 0; j < n; ++j ) f += x(j);
    f = exp(cos((k+1)*f))-1;
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
  	real_type acc = 0;
    for ( int_type j = 0; j < n; ++j ) acc += x(j);
    for ( int_type j = 0; j < n; ++j ) f(j) = exp(cos((j+1)*acc));
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
  	real_type acc = 0;
    for ( int_type j = 0; j < n; ++j ) acc += x(j);

    int_type kk = 0;
    for ( int_type i = 0; i < n; ++i ) {
      real_type cc = cos((i+1)*acc);
      real_type ss = sin((i+1)*acc);
      for ( int_type j = 0; j < n; ++j )
        jac(kk++) = -exp(cc)*ss*(i+1);
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k ) x(k) = 0;
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

  string note() const { return "This problem has no solution"; }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo211 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo211( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.211",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = k+1;
    for ( int_type j = 0; j < n; ++j ) f += power3(x(j));
    return f/(2*n);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
  	real_type acc = 0;
    for ( int_type j = 0; j < n; ++j ) acc += power3(x(j));
    for ( int_type j = 0; j < n; ++j ) f(j) = (acc+j+1)/(2*n);
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
  	real_type tmp = 1.5/n;
    for ( int_type i = 0; i < n; ++i )
      for ( int_type j = 0; j < n; ++j )
        jac(kk++) = tmp*power2(x(j));
  }

  virtual
  int_type
  numExactSolution() const override
  { return 0; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override
  { }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.setZero();
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

  string note() const { return "This problem has no solution"; }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo212 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo212( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.212",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f = 0;
    if ( k == 0 ) f = x(0);
    else          f = cos(x(k-1))+x(k)-1;
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = x(0);
    for ( int_type k = 1; k < n; ++k )
      f(k) = cos(x(k-1))+x(k)-1;
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 2*n-1;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    for ( int_type k = 1; k < n; ++k ) {
      SETIJ(k,k-1);
      SETIJ(k,k);
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    jac(kk++) = 1;
    for ( int_type k = 1; k < n; ++k ) {
      jac(kk++) = -sin(x(k-1));
      jac(kk++) = 1;
    }
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }

  virtual
  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x.setZero();
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(0.5);
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo213 : public nonlinearSystem {
  real_type h2;
public:
  
  RooseKullaLombMeressoo213( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.213",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); h2 = 1.0/power2(n+1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type xp = k < n-1 ? x(k+1) : 1;
  	real_type xm = k > 0   ? x(k-1) : 0;
  	real_type xc = x(k);
    return 2*xc-xm-xp+h2*(xc+sin(xc));
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k )
      f(k) = evalFk( x, k );
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 0; k < n; ++k ) {
      if ( k > 0 ) { SETIJ(k,k-1); }
      SETIJ(k,k);
      if ( k < n-1 ) { SETIJ(k,k+1); }
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n; ++k ) {
      if ( k > 0 ) jac(kk++) = -1;
      jac(kk++) = 2+h2*(1+cos(x(k)));
      if ( k < n-1 ) jac(kk++) = -1;
    }
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 1;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = 0.25493102645914878;
      x(1) = 0.566207574177154838;
      break;
    case 5:
      x(0) = 0.123717884647979601;
      x(1) = 0.254300224963985633;
      x(2) = 0.398934465982092534;
      x(3) = 0.565440128266567754;
      x(4) = 0.762535446654485805;
      break;
    case 10:
      x(0) = 0.0670031343864649215;
      x(1) = 0.135113344358246179;
      x(2) = 0.20545343835958349;
      x(3) = 0.279177536667768345;
      x(4) = 0.35748628505027491;
      x(5) = 0.441641370870779482;
      x(6) = 0.532978814497350184;
      x(7) = 0.632920225489453259;
      x(8) = 0.74298082731276216;
      x(9) = 0.864772561262492578;
      break;
    case 20:
      x(0) = 0.0350196898973525322;
      x(1) = 0.0701981830199509271;
      x(2) = 0.105694904547656857;
      x(3) = 0.141670522051372094;
      x(4) = 0.178287562727141746;
      x(5) = 0.21571102540715037;
      x(6) = 0.254108984809966743;
      x(7) = 0.293653184781339449;
      x(8) = 0.334519616353437577;
      x(9) = 0.376889075292977049;
      x(10) = 0.42094769239294999;
      x(11) = 0.466887428063573084;
      x(12) = 0.514906520772035803;
      x(13) = 0.565209876549352397;
      x(14) = 0.618009384118297755;
      x(15) = 0.673524137209132956;
      x(16) = 0.731980542357675512;
      x(17) = 0.793612287002704475;
      x(18) = 0.858660139154802993;
      x(19) = 0.927371546514711986;
      break;
    case 30:
      x(0) = 0.0237123159942052714;
      x(1) = 0.0474739789290085645;
      x(2) = 0.0713344245144125061;
      x(3) = 0.095343265901786059;
      x(4) = 0.119550382155985915;
      x(5) = 0.144006006415265148;
      x(6) = 0.168760813611468746;
      x(7) = 0.193866007602538554;
      x(8) = 0.21937340754326351;
      x(9) = 0.245335533288195989;
      x(10) = 0.271805689582330812;
      x(11) = 0.298838048750088958;
      x(12) = 0.32648773154088373;
      x(13) = 0.354810885729579639;
      x(14) = 0.383864762001958348;
      x(15) = 0.413707786578368653;
      x(16) = 0.444399629942578145;
      x(17) = 0.47600127094703848;
      x(18) = 0.50857505546002324;
      x(19) = 0.542184748604255695;
      x(20) = 0.576895579510812273;
      x(21) = 0.612774277376699539;
      x(22) = 0.649889097470416455;
      x(23) = 0.688309835578470475;
      x(24) = 0.728107829229357573;
      x(25) = 0.769355943873014314;
      x(26) = 0.812128542037363621;
      x(27) = 0.856501433334843409;
      x(28) = 0.902551803057848212;
      x(29) = 0.950358116991871227;
      break;
    }
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(1);
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo214 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo214( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.214",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type f    = 0;
    int_type  imin = std::max(0,k-5);
    int_type  imax = std::min(n-1,k+1);
    for ( int_type i = imin; i <= imax; ++i )
      if ( i != k )
        f += x(i)+power2(x(i));
    return x(k)*(2+5*power2(x(k))) + 1 - f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type j = 0; j < n; ++j )
      f(j) = evalFk(x,j);
  }

  virtual
  int_type
  jacobianNnz() const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n; ++k ) {
      int_type imin = std::max(0,k-5);
      int_type imax = std::min(n-1,k+1);
      kk += imax - imin + 1;
    }
    return kk;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 0; k < n; ++k ) {
      int_type imin = std::max(0,k-5);
      int_type imax = std::min(n-1,k+1);
      for ( int_type i = imin; i <= imax; ++i ) {
        SETIJ(k,i);
      }
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n; ++k ) {
      int_type imin = std::max(0,k-5);
      int_type imax = std::min(n-1,k+1);
      for ( int_type i = imin; i <= imax; ++i ) {
        if ( i == k ) jac(kk) = 2+15*power2(x(k));
        else          jac(kk) = -(1+2*x(i));
        ++kk;
      }
    }
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 1;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = -0.427304623558166286;
      x(1) = -0.427304623558166286;
      break;
    case 5:
      x(0) = -0.428302864642700787;
      x(1) = -0.476596531501095377;
      x(2) = -0.519637722100754651;
      x(3) = -0.558861956527025194;
      x(4) = -0.558861956527025194;
      break;
    case 10:
      x(0) = -0.428302863587250282;
      x(1) = -0.476596424356290238;
      x(2) = -0.519652463646861684;
      x(3) = -0.558099324832180832;
      x(4) = -0.592506156829457287;
      x(5) = -0.624503682199467947;
      x(6) = -0.623239471440591108;
      x(7) = -0.621393841796573532;
      x(8) = -0.620453596659087392;
      x(9) = -0.586469270720434976;
      break;
    case 20:
      x(0) = -0.428302863587250282;
      x(1) = -0.476596424356293569;
      x(2) = -0.519652463646401386;
      x(3) = -0.55809932485615199;
      x(4) = -0.592506155965082826;
      x(5) = -0.624503707410516529;
      x(6) = -0.623238669132451184;
      x(7) = -0.621419676713647839;
      x(8) = -0.61961584283347626;
      x(9) = -0.618226017919857429;
      x(10) = -0.617518024841495761;
      x(11) = -0.617731830318644759;
      x(12) = -0.617900316253351289;
      x(13) = -0.618007798540867848;
      x(14) = -0.618057061755049264;
      x(15) = -0.618062699716298014;
      x(16) = -0.618047199350808651;
      x(17) = -0.618011195738616514;
      x(18) = -0.618872079495047522;
      x(19) = -0.586276945400115101;
      break;
    case 30:
      x(0) = -0.428302863587250282;
      x(1) = -0.476596424356293569;
      x(2) = -0.519652463646401386;
      x(3) = -0.55809932485615199;
      x(4) = -0.592506155965082826;
      x(5) = -0.624503707410516529;
      x(6) = -0.623238669132451184;
      x(7) = -0.621419676713647839;
      x(8) = -0.61961584283347626;
      x(9) = -0.618226017919857429;
      x(10) = -0.617518024841495206;
      x(11) = -0.617731830318665742;
      x(12) = -0.617900316252663617;
      x(13) = -0.61800779856335919;
      x(14) = -0.618057061019478993;
      x(15) = -0.618062723774471579;
      x(16) = -0.61804641236762925;
      x(17) = -0.618036943255954929;
      x(18) = -0.618032796823900332;
      x(19) = -0.618032010907616169;
      x(20) = -0.618032748437421176;
      x(21) = -0.61803365220978157;
      x(22) = -0.618034039196207474;
      x(23) = -0.618034129052205672;
      x(24) = -0.618034091025163379;
      x(25) = -0.618034003909173957;
      x(26) = -0.618034776213912562;
      x(27) = -0.618008230615912701;
      x(28) = -0.618873272626757731;
      x(29) = -0.58627911806458255;
      break;
    }
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(-1);
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo215 : public nonlinearSystem {
  mutable dvec_t grad_dS;
  mutable dvec_t grad_S;

public:
  
  RooseKullaLombMeressoo215( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.215",RKM_BIBTEX,neq)
  { checkMinEquations(n,1); grad_dS.resize(n); grad_S.resize(n); }

  real_type
  g( int_type k ) const
  { return (k+1.0)/29.0; }

  real_type
  S( dvec_t const & x, int_type k ) const {
    real_type sum = 0;
    real_type gkj = 1;
    real_type gk  = g(k);
    for ( int_type j = 0; j < n; ++j, gkj *= gk ) sum += gkj * x(j);
    return sum;
  }

  void
  S_1( dvec_t const & x, int_type k, dvec_t & grad ) const {
    real_type gkj = 1;
    real_type gk  = g(k);
    for ( int_type j = 0; j < n; ++j, gkj *= gk ) grad[j] = gkj;
  }

  real_type
  dS( dvec_t const & x, int_type k ) const {
    real_type sum = 0;
    real_type gkj = 1;
    real_type gk  = g(k);
    for ( int_type j = 1; j < n; ++j, gkj *= gk ) sum += j * gkj * x(j);
    return sum;
  }

  void
  dS_1( dvec_t const & x, int_type k, dvec_t & grad ) const {
    real_type gkj = 1;
    real_type gk  = g(k);
    grad[0] = 0;
    for ( int_type j = 1; j < n; ++j, gkj *= gk ) grad[j] = j * gkj;
  }
  
  real_type
  powergk( int_type k, int_type i ) const {
    real_type gk  = g(k);
    if ( i == 0 ) return 1/gk;
    real_type res = 1;
    for ( int_type j = 1; j < i; ++j ) res *= gk;
    return res;
  }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
  	real_type f = 0;
    for ( int_type k = 0; k < 29; ++k ) {
      real_type Sk    = S(x,k);
      real_type dSk   = dS(x,k);
      real_type gk    = g(k);
      real_type powgk = powergk(k,i);
      f += powgk*(i-2*gk*Sk)*(dSk-Sk*Sk-1);
    }
    if ( i == 0 ) f += x(0)*(1-2*(x(1)-x(0)*x(0)-1));
    if ( i == 1 ) f += x(1)-x(0)*x(0)-1;
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type j = 0; j < n; ++j )
      f(j) = evalFk(x,j);
  }

  virtual
  int_type
  jacobianNnz() const override
  { return n*n; }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0; // fortran addressing
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    for ( int_type i = 0; i < n; ++i ) {
      for ( int_type k = 0; k < 29; ++k ) {
        real_type gk    = g(k);
        real_type powgk = powergk(k,i);
        real_type Sk    = S(x,k);
        real_type dSk   = dS(x,k);
        dS_1( x, k, grad_dS );
        S_1( x, k, grad_S );
        real_type A = i-2*gk*Sk;
        real_type B = dSk-Sk*Sk-1;
        for ( int_type j = 0; j < n; ++j ) {
          real_type A_1 = -2*gk*grad_S[j];
          real_type B_1 = grad_dS[j]-2*Sk*grad_S[j];
          jac[caddr(i,j)] += powgk*(A_1*B+A*B_1);
        }
      }
    }
    jac[caddr(0,0)] += 6*x(0)*x(0)-2*x(1)+3;
    jac[caddr(0,1)] -= 2*x(0);
    jac[caddr(1,0)] -= 2*x(0);
    jac[caddr(1,1)] += 1;
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 6:
    case 10:
      return 1;
    default:
    case 20:
    case 30:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = -0.501367007521962726;
      x(1) = 1.0736498384595432;
      break;
    case 5:
      x(0) = -0.0714051701051878485;
      x(1) = 0.970024704720390041;
      x(2) = 0.2661688465197945;
      x(3) = -0.539893187056864843;
      x(4) = 0.710711221847347585;
      break;
    case 6:
      x(0) = -0.0157250864014588515;
      x(1) = 1.01243486936910987;
      x(2) = -0.232991625956739556;
      x(3) = 1.26043008779961352;
      x(4) = -1.51372892272228565;
      x(5) = 0.992996432431136333;
      break;
    case 10:
      x(0) = -1.22248986826704962e-06;
      x(1) = 1.00004450564633962;
      x(2) = -0.00508490059336646257;
      x(3) = 0.417429340147828398;
      x(4) = -0.588443484856691734;
      x(5) = 2.30197280176861696;
      x(6) = -4.56245409953258463;
      x(7) = 5.59197371897036266;
      x(8) = -3.64040714755196371;
      x(9) = 1.04237104657441626;
      break;
    case 20:
      break;
    case 30:
      break;
    }
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.setZero();
    //getExactSolution( x, 0 );
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( x(i) > -50 && x(i) < 50, "Bad range" );
  }

  string note() const { return "Watson Function"; }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo216 : public nonlinearSystem {
  mutable std::vector<std::vector<real_type> > m_y, m_y_D;
public:
  
  RooseKullaLombMeressoo216( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.216",RKM_BIBTEX,neq) {
    NONLIN_ASSERT(
      n > 0 && n != 8 && n < 10,
      "RooseKullaLombMeressoo216, neq=" << n << " must be [1..7] or 9"
    );
    m_y.resize(n);
    m_y_D.resize(n);
    for ( int_type i = 0; i < n; ++i ) {
      m_y[i].resize(n);
      m_y_D[i].resize(n);
    }
  }

  void
  evalY( real_type s, std::vector<real_type> & y ) const {
    y[0] = 2*s-1;
    y[1] = 8*(s-1)*s+1;
    for ( int_type k = 2; k < n; ++k )
      y[k] = 2*y[0]*y[k-1]-y[k-2];
  }

  void
  evalY_D(
    real_type s,
    std::vector<real_type> & y,
    std::vector<real_type> & y_D
  ) const {
    evalY( s, y );
    y_D[0] = 2;
    y_D[1] = 16*s-8;
    for ( int_type k = 2; k < n; ++k )
      y_D[k] = 2*(y[0]*y_D[k-1]+y_D[0]*y[k-1])-y_D[k-2];
  }

  void
  evalY( dvec_t const & x ) const {
    for ( int_type k = 0; k < n; ++k )
      evalY( x(k), m_y[k] );
  }

  void
  evalY_D( dvec_t const & x ) const {
    for ( int_type k = 0; k < n; ++k )
      evalY_D( x(k), m_y[k], m_y_D[k] );
  }

  real_type
  Y( int_type xj, int_type k ) const
  { return m_y[xj][k]; }

  real_type
  Y_D( int_type xj, int_type k ) const
  { return m_y_D[xj][k]; }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type f = 0;
    evalY( x );
    for ( int_type j = 0; j < n; ++j ) f += Y(j,k);
    f /= n;
    if ( (k % 2) == 1 ) f += 1.0/(power2(k+1)-1.0);
    return f;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    evalY( x );
    for ( int_type k = 0; k < n; ++k ) {
      f(k) = 0;
      for ( int_type j = 0; j < n; ++j ) f(k) += Y(j,k);
      f(k) /= n;
      if ( (k % 2) == 1 ) f(k) += 1.0/(power2(k+1)-1.0);
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
    evalY_D( x );
    int_type kk = 0;
    for ( int_type k = 0; k < n; ++k )
      for ( int_type j = 0; j < n; ++j )
        jac(kk++) = Y_D(j,k)/n;
  }

  virtual
  int_type
  numExactSolution() const override
  { return 1; }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = 0.211324865405187134;
      x(1) = 0.788675134594812866;
      break;
    case 3:
      x(0) = 0.146446609406726241;
      x(1) = 0.5;
      x(2) = 0.853553390593273731;
      break;
    case 4:
      x(0) = 0.10267276385411693;
      x(1) = 0.406203762957460079;
      x(2) = 0.593796237042539921;
      x(3) = 0.897327236145883056;
      break;
    case 5:
      x(0) = 0.0837512564995090553;
      x(1) = 0.312729295223209469;
      x(2) = 0.5;
      x(3) = 0.687270704776790531;
      x(4) = 0.916248743500490903;
      break;
    case 6:
      x(0) = 0.0668765909460896923;
      x(1) = 0.288740673119444125;
      x(2) = 0.366682299241647691;
      x(3) = 0.633317700758352253;
      x(4) = 0.71125932688055582;
      x(5) = 0.933123409053910335;
      break;
    case 7:
      x(0) = 0.0580691496209754937;
      x(1) = 0.23517161235742165;
      x(2) = 0.338044094740046208;
      x(3) = 0.5;
      x(4) = 0.661955905259953847;
      x(5) = 0.764828387642578433;
      x(6) = 0.941930850379024465;
      break;
    case 9:
      x(0) = 0.0442053461357827734;
      x(1) = 0.199490672309881045;
      x(2) = 0.235619108471059907;
      x(3) = 0.416046907892598183;
      x(4) = 0.499999999999999833;
      x(5) = 0.583953092107402094;
      x(6) = 0.764380891528939954;
      x(7) = 0.800509327690119066;
      x(8) = 0.955794653864217247;
      break;
    }
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    real_type bf = 1.0/(n+1);
    for ( int_type k = 0; k < n; ++k ) x(k) = bf*(k+1);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override {
    for ( int_type i = 0; i < n; ++i )
      NONLIN_ASSERT( x(i) > 0, "Bad range" );
  }

  string note() const { return "Each permutation of x is a solution"; }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo217 : public nonlinearSystem {
public:
  
  RooseKullaLombMeressoo217( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.217",RKM_BIBTEX,neq)
  { checkMinEquations(n,2); }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type xc = x(k);
  	real_type xp = k < n-1 ? x(k+1) : 20;
  	real_type xm = k > 0   ? x(k-1) : 0;
    return 3*xc*(xm+xp-2*xc)+power2(xp-xm)/4;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k )
      f(k) = evalFk( x, k );
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 0; k < n; ++k ) {
      SETIJ(k,k);
      if ( k > 0   ) { SETIJ(k,k-1); }
      if ( k < n-1 ) { SETIJ(k,k+1); }
    }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 0; k < n; ++k ) {
    	real_type xc = x(k);
    	real_type xp = k < n-1 ? x(k+1) : 20;
    	real_type xm = k > 0   ? x(k-1) : 0;
      jac(kk++) = 3*(xp+xm-4*xc);
      if ( k > 0   ) jac(kk++) = 3*xc-(xp-xm)/2;
      if ( k < n-1 ) jac(kk++) = 3*xc+(xp-xm)/2;
    }
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 1;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = 8.33828484333201025;
      x(1) = 14.5583676083251294;
      break;
    case 5:
      x(0) = 4.88928869829531543;
      x(1) = 8.53653521682389993;
      x(2) = 11.7273247200356803;
      x(3) = 14.6523320192791218;
      x(4) = 17.394667732103354;
      break;
    case 10:
      x(0) = 3.08315248959636978;
      x(1) = 5.38308155447113368;
      x(2) = 7.39517190291697979;
      x(3) = 9.23966178544206862;
      x(4) = 10.9689601971417012;
      x(5) = 12.6118651601462197;
      x(6) = 14.1863707080987851;
      x(7) = 15.7046865038081975;
      x(8) = 17.1755885168751945;
      x(9) = 18.6056591191924703;
      break;
    case 20:
      x(0) = 1.89123927553494386;
      x(1) = 3.30204078247077248;
      x(2) = 4.53627888965039894;
      x(3) = 5.66770904788268393;
      x(4) = 6.72847950486216018;
      x(5) = 7.73625528062737722;
      x(6) = 8.70207411115346829;
      x(7) = 9.63342553642360322;
      x(8) = 10.53569283165554;
      x(9) = 11.4129136953710457;
      x(10) = 12.2682175592459135;
      x(11) = 13.1040939164598367;
      x(12) = 13.9225656143308374;
      x(13) = 14.7253054956124636;
      x(14) = 15.5137176336911349;
      x(15) = 16.2889955537418949;
      x(16) = 17.0521649909861708;
      x(17) = 17.8041159610068647;
      x(18) = 18.5456272591407938;
      x(19) = 19.2773854806811933;
      break;
    case 30:
      x(0) = 1.41029144531941619;
      x(1) = 2.46232189012527636;
      x(2) = 3.38268953823807861;
      x(3) = 4.22639360771921346;
      x(4) = 5.01740695028136763;
      x(5) = 5.76890231829069311;
      x(6) = 6.48911051829969399;
      x(7) = 7.18361647777063173;
      x(8) = 7.85643448885871631;
      x(9) = 8.51057545121989989;
      x(10) = 9.14837297265207816;
      x(11) = 9.77168346073938032;
      x(12) = 10.3820153466489788;
      x(13) = 10.9806160641959405;
      x(14) = 11.5685326266858191;
      x(15) = 12.1466550422555599;
      x(16) = 12.7157481985778418;
      x(17) = 13.2764757775986784;
      x(18) = 13.8294185246607579;
      x(19) = 14.3750884318112888;
      x(20) = 14.9139399076884889;
      x(21) = 15.4463786872421789;
      x(22) = 15.9727690205710928;
      x(23) = 16.4934395336305428;
      x(24) = 17.0086880512958594;
      x(25) = 17.5187856006541267;
      x(26) = 18.0239797600299596;
      x(27) = 18.5244974809442624;
      x(28) = 19.0205474818159175;
      x(29) = 19.5123222909239296;
      break;
    }
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x.fill(10);
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

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo218 : public nonlinearSystem {
  real_type h, h2;
public:
  
  RooseKullaLombMeressoo218( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.218",RKM_BIBTEX,neq)
  { checkMinEquations(n,2); h = 1.0/(n+1); h2 = h*h; }
  
  real_type
  t( int_type j ) const
  { return (j+1)*h; }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
  	real_type xc = x(k);
  	real_type xp = k < n-1 ? x(k+1) : 0;
  	real_type xm = k > 0   ? x(k-1) : 0;
    return 2*xc-(xm+xp)+(0.5*h2)*power3(xc+t(k)+1);
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k )
      f(k) = evalFk( x, k );
  }

  virtual
  int_type
  jacobianNnz() const override {
    return 3*n-2;
  }

  virtual
  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( int_type k = 1; k < n;   ++k ) { SETIJ(k,k-1); }
    for ( int_type k = 0; k < n-1; ++k ) { SETIJ(k,k+1); }
    for ( int_type k = 0; k < n;   ++k ) { SETIJ(k,k); }
    #undef SETIJ
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    int_type kk = 0;
    for ( int_type k = 1; k < n;   ++k ) jac(kk++) = -1;
    for ( int_type k = 0; k < n-1; ++k ) jac(kk++) = -1;
    for ( int_type k = 0; k < n;   ++k ) jac(kk++) = 2+1.5*h2*power2(x(k)+t(k)+1);
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 1;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = -0.128246763033731614;
      x(1) = -0.159267567244640834;
      break;
    case 5:
      x(0) = -0.0750221292923204525;
      x(1) = -0.131976210352190593;
      x(2) = -0.164848771909337249;
      x(3) = -0.164664680215800663;
      x(4) = -0.117417651684193575;
      break;
    case 10:
      x(0) = -0.0431649825187648689;
      x(1) = -0.0815771565353868716;
      x(2) = -0.114485714380529277;
      x(3) = -0.14097357686259665;
      x(4) = -0.159908696181983112;
      x(5) = -0.169877202312774922;
      x(6) = -0.169089983781208347;
      x(7) = -0.155249535221831853;
      x(8) = -0.125355891678934989;
      x(9) = -0.0754165336858920871;
      break;
    case 20:
      x(0) = -0.0232104399966753146;
      x(1) = -0.0452020277196446968;
      x(2) = -0.0658809801847738824;
      x(3) = -0.0851436508063508207;
      x(4) = -0.102875197721159828;
      x(5) = -0.118948030435962207;
      x(6) = -0.133219990530385718;
      x(7) = -0.145532211750480533;
      x(8) = -0.155706591602235406;
      x(9) = -0.16354278961733093;
      x(10) = -0.16881464562153628;
      x(11) = -0.171265882959322213;
      x(12) = -0.170604924470161234;
      x(13) = -0.166498599947109888;
      x(14) = -0.158564458442799261;
      x(15) = -0.146361310877112738;
      x(16) = -0.129377508964773374;
      x(17) = -0.107016302444755948;
      x(18) = -0.0785773886624082651;
      x(19) = -0.0432334478445528594;
      break;
    case 30:
      x(0) = -0.0158588747608703062;
      x(1) = -0.0311714390223493919;
      x(2) = -0.0459099102817520852;
      x(3) = -0.060044590713602998;
      x(4) = -0.073543699225747064;
      x(5) = -0.0863731855306687224;
      x(6) = -0.0984965239444887258;
      x(7) = -0.10987448428747057;
      x(8) = -0.120464876863777437;
      x(9) = -0.130222268033639288;
      x(10) = -0.139097662344602108;
      x(11) = -0.147038146543775367;
      x(12) = -0.153986490029936474;
      x(13) = -0.159880695398378264;
      x(14) = -0.164653491652127226;
      x(15) = -0.168231761362994642;
      x(16) = -0.170535891518048321;
      x(17) = -0.17147903592302513;
      x(18) = -0.17096627478051718;
      x(19) = -0.168893654324848153;
      x(20) = -0.165147086059982073;
      x(21) = -0.159601081061931882;
      x(22) = -0.15211728978118419;
      x(23) = -0.142542811566464672;
      x(24) = -0.130708230408252885;
      x(25) = -0.116425323750638757;
      x(26) = -0.0994843790941261352;
      x(27) = -0.0796510377833258565;
      x(28) = -0.0566625658741519017;
      x(29) = -0.0302234270054019573;
      break;
    }
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k ) x(k) = t(k)*(t(k)-1);
  }

  virtual
  int_type
  numInitialPoint() const override
  { return 1; }

  virtual
  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

  string note() const { return "For any n all solutions are in the interval [-0.5,0]"; }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class RooseKullaLombMeressoo219 : public nonlinearSystem {
  real_type h;
public:
  
  RooseKullaLombMeressoo219( int_type neq )
  : nonlinearSystem("Roose Kulla Lomb Meressoo N.219",RKM_BIBTEX,neq)
  { checkMinEquations(n,2); h = 1.0/(n+1); }
  
  real_type
  t( int_type j ) const
  { return (j+1)*h; }

  virtual
  real_type
  evalFk( dvec_t const & x, int_type i ) const override {
  	real_type sum1 = 0;
  	real_type sum2 = 0;
    for ( int_type j = 0;   j <= i; ++j ) sum1 += t(j)*power3(x(j)+t(j)+1);
    for ( int_type j = i+1; j < n;  ++j ) sum2 += (1-t(j))*power3(x(j)+t(j)+1);
    return x(i)+0.5*h*( (1-t(i))*sum1 + t(i)*sum2 );
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( int_type k = 0; k < n; ++k )
      f(k) = evalFk(x,k);
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
    
      real_type bf = 1.5*h*(1-t(i));
      for ( int_type j = 0; j <= i; ++j )
        jac(kk++) = bf*t(j)*power2(x(j)+t(j)+1);

      jac[kk-1] += 1; // somma su j(i,i)

      bf = 1.5*h*t(i);
      for ( int_type j = i+1; j < n; ++j )
        jac(kk++) = bf*(1-t(j))*power2(x(j)+t(j)+1);
    }
  }

  virtual
  int_type
  numExactSolution() const override {
    switch ( n ) {
    case 2:
    case 5:
    case 10:
    case 20:
    case 30:
      return 1;
    default:
      return 0;
    }
  }
 
  virtual
  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
    switch ( n ) {
    case 2:
      x(0) = -0.128246763033731614;
      x(1) = -0.159267567244640834;
      break;
    case 5:
      x(0) = -0.0750221292923204663;
      x(1) = -0.13197621035219062;
      x(2) = -0.164848771909337277;
      x(3) = -0.164664680215800718;
      x(4) = -0.117417651684193616;
      break;
    case 10:
      x(0) = -0.043164982518764862;
      x(1) = -0.0815771565353868855;
      x(2) = -0.114485714380529277;
      x(3) = -0.14097357686259665;
      x(4) = -0.159908696181983112;
      x(5) = -0.169877202312774894;
      x(6) = -0.169089983781208347;
      x(7) = -0.155249535221831825;
      x(8) = -0.125355891678934933;
      x(9) = -0.0754165336858920177;
      break;
    case 20:
      x(0) = -0.0232104399966752729;
      x(1) = -0.0452020277196446066;
      x(2) = -0.0658809801847737575;
      x(3) = -0.0851436508063506542;
      x(4) = -0.102875197721159606;
      x(5) = -0.118948030435961943;
      x(6) = -0.133219990530385413;
      x(7) = -0.145532211750480228;
      x(8) = -0.1557065916022351;
      x(9) = -0.163542789617330597;
      x(10) = -0.168814645621535919;
      x(11) = -0.171265882959321825;
      x(12) = -0.170604924470160901;
      x(13) = -0.166498599947109527;
      x(14) = -0.158564458442798956;
      x(15) = -0.146361310877112488;
      x(16) = -0.129377508964773152;
      x(17) = -0.107016302444755837;
      x(18) = -0.0785773886624082513;
      x(19) = -0.0432334478445528664;
      break;
    case 30:
      x(0) = -0.0158588747608703166;
      x(1) = -0.0311714390223494162;
      x(2) = -0.0459099102817521199;
      x(3) = -0.0600445907136030396;
      x(4) = -0.0735436992257471334;
      x(5) = -0.0863731855306687918;
      x(6) = -0.0984965239444888091;
      x(7) = -0.109874484287470653;
      x(8) = -0.120464876863777492;
      x(9) = -0.130222268033639343;
      x(10) = -0.139097662344602191;
      x(11) = -0.147038146543775394;
      x(12) = -0.153986490029936474;
      x(13) = -0.159880695398378181;
      x(14) = -0.164653491652127115;
      x(15) = -0.168231761362994447;
      x(16) = -0.170535891518048127;
      x(17) = -0.171479035923024881;
      x(18) = -0.170966274780516903;
      x(19) = -0.168893654324847792;
      x(20) = -0.165147086059981713;
      x(21) = -0.159601081061931466;
      x(22) = -0.152117289781183801;
      x(23) = -0.142542811566464256;
      x(24) = -0.130708230408252524;
      x(25) = -0.116425323750638438;
      x(26) = -0.0994843790941258299;
      x(27) = -0.0796510377833256483;
      x(28) = -0.0566625658741517907;
      x(29) = -0.0302234270054018567;
      break;
    }
  }
  
  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    for ( int_type k = 0; k < n; ++k ) x(k) = t(k)*(t(k)-1);
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
