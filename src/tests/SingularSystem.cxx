/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

#define SINGULAR_SYSTEM_BIBTEX \
"@article{Waziri:2011,\n" \
"  author  = {Waziri Yusuf, Mohammed and June, Leong Wah and Hassan, Malik Abu},\n" \
"  title   = {Jacobian-free diagonal {N}ewton's method for solving nonlinear\n" \
"             systems with singular {J}acobian},\n" \
"  joirnal = {Malaysian Journal of Mathematical Sciences},\n" \
"  volume  = {5},\n" \
"  year    = {2011},\n" \
"  number  = {2},\n" \
"  pages   = {241--255}\n" \
"}\n"

/*
  ￼Modified Newton’s method for systems of nonlinear equations with singular Jacobian
*/

// Problem 1. (Jose et al. (2009))

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemA : public nonlinearSystem {
public:

  SingularSystemA()
  : nonlinearSystem("Singular System A",SINGULAR_SYSTEM_BIBTEX,2)
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    if ( i == 0 ) return power2(x(0)-1)*(x(0)-x(1));
    else          return power5(x(1)-2)*cos(2*x(0)/x(1));
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(x(0)-1)*(x(0)-x(1));
    f(1) = power5(x(1)-2)*cos(2*x(0)/x(1));
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type t = 2*x(0)/x(1);
    real_type S = sin(t);
    real_type C = cos(t);
    jac(0) = (x(0)-1)*(3*x(0)-2*x(1)-1);
    jac(1) = -power2(x(0)-1);
    jac(2) = -2*power5(x(1)-2)*S/x(1);
    jac(3) = power4(x(1)-2)*(5*C+2*(x(1)-2)*x(0)*S/(x(1)*x(1)));
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 2;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 1.5; x(1) = 2.5; break;
      case 1: x(0) = 2;   x(1) = 5;   break;
      case 2: x(0) = 0;   x(1) = 3;   break;
      case 3: x(0) = 0.5; x(1) = 2;   break;
    }
  }

  integer
  numInitialPoint() const override
  { return 4; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemB : public nonlinearSystem {
public:

  SingularSystemB()
  : nonlinearSystem("Singular System B",SINGULAR_SYSTEM_BIBTEX,3)
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & X, dvec_t & f ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);
    f(0) = power4(x-1)*exp(y);
    f(1) = power5(y-2)*(x*y-1);
    f(2) = power6(z+4);
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
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);
    jac(0) = 4*power3(x-1)*exp(y);
    jac(1) = power4(x-1)*exp(y);
    jac(2) = 0;
    jac(3) = power5(y-2)*y;
    jac(4) = power4(y-2)*(6*x*y-2*x-5);
    jac(5) = 0;
    jac(6) = 0;
    jac(7) = 0;
    jac(8) = 6*power5(z+4);
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 2;
    x(2) = -4;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 2; x(1) = 1; x(2) = -2; break;
      case 1: x(0) = 4; x(1) = 5; x(2) =  6; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemC : public nonlinearSystem {
public:

  SingularSystemC()
  : nonlinearSystem("Singular System C",SINGULAR_SYSTEM_BIBTEX,2)
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & X, dvec_t & f ) const override {
    real_type x = X(0);
    real_type y = X(1);
    f(0) = exp(x)-y-1;
    f(1) = x-y;
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    real_type x = X(0);
    jac(0) = exp(x);
    jac(1) = -1;
    jac(2) = 1;
    jac(3) = -1;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0.7; x(1) = 0.7; break;
      case 1: x(0) = 2;   x(1) = 1;   break;
    }
  }

  integer
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemD : public nonlinearSystem {
public:

  SingularSystemD()
  : nonlinearSystem("Singular System D",SINGULAR_SYSTEM_BIBTEX,2)
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & X, dvec_t & f ) const override {
    real_type x = X(0);
    real_type y = X(1);
    f(0) = power4(6*x-y);
    f(1) = cos(x)-1+y;
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    real_type x = X(0);
    real_type y = X(1);
    jac(0) = 24*power3(6*x-y);
    jac(1) = -4*power3(6*x-y);
    jac(2) = -sin(x);
    jac(3) = 1;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = -0.5; x(1) = -0.5; break;
      case 1: x(0) = 4;    x(1) = 1;    break;
    }
  }

  integer
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemE : public nonlinearSystem {
public:

  SingularSystemE()
  : nonlinearSystem("Singular System E",SINGULAR_SYSTEM_BIBTEX,3)
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & X, dvec_t & f ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);
    f(0) = 3*x-cos(y*z)-0.5;
    f(1) = x*x-635*y*y-0.25;
    f(2) = exp(-x*y)+20*z+(10*m_pi-3)/3;
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
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);
    jac(0) = 3;
    jac(1) = z*sin(y*z);
    jac(2) = y*sin(y*z);
    jac(3) = 2*x;
    jac(4) = -1270*y;
    jac(5) = 0;
    jac(6) = -y*exp(-x*y);
    jac(7) = -x*exp(-x*y);
    jac(8) = 20;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0.5;
    x(1) = 0;
    x(2) = -m_pi/6;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0.2; x(1) = 0.2; x(2) = -0.2; break;
      case 1: x(0) = 1;   x(1) = 2;   x(2) =  1;   break;
    }
  }

  integer
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemF : public nonlinearSystem {
public:

  SingularSystemF()
  : nonlinearSystem("Singular System F",SINGULAR_SYSTEM_BIBTEX,3)
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & X, dvec_t & f ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);
    f(0) = exp(x*x)-8*x*sin(y);
    f(1) = x+y-1;
    f(2) = power3(z-1);
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
  jacobian( dvec_t const & X, dvec_t & jac ) const override {
    real_type x = X(0);
    real_type y = X(1);
    real_type z = X(2);
    jac(0) = 2*x*exp(x*x)-8*sin(y);
    jac(1) = -8*x*cos(y);
    jac(2) = 0;
    jac(3) = 1;
    jac(4) = 1;
    jac(5) = 0;
    jac(6) = 0;
    jac(7) = 0;
    jac(8) = 3*power2(z-1);
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0.175599;
    x(1) = 0.824401;
    x(2) = 1;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0; x(1) = 1; x(2) = 2; break;
      case 1: x(0) = 1; x(1) = 1; x(2) = 3; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP2 : public nonlinearSystem {
public:

  SingularSystemP2()
  : nonlinearSystem(
      "Singular System Problem 2 (Ishihara, K. 2001)",
      SINGULAR_SYSTEM_BIBTEX,
      3
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 4*x(0)-2*x(1)+power2(x(0))-3;
    f(1) = -x(0)+4*x(1)-x(2)+power2(x(1))-3;
    f(2) = -2*x(1)+4*x(2)+power2(x(2))-3;
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
    jac(0) = 4+2*x(0);
    jac(1) = -2;
    jac(2) = 0;

    jac(3) = -1;
    jac(4) = 4+2*x(1);
    jac(5) = -1;

    jac(6) = 0;
    jac(7) = -2;
    jac(8) = 4+2*x(2);
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 1;
    x(2) = 1;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = -1.5; x(1) = 0; x(2) = -1.5; break;
      case 1: x(0) =    4; x(1) = 0; x(2) = 4;    break;
      case 2: x(0) =   -1; x(1) = 5; x(2) = -1;   break;
      case 3: x(0) =    4; x(1) = 4; x(2) = 4;    break;
      case 4: x(0) =  -10; x(1) = 0; x(2) = 10;   break;
    }
  }

  integer
  numInitialPoint() const override
  { return 5; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP3 : public nonlinearSystem {
public:

  SingularSystemP3()
  : nonlinearSystem(
      "Singular System Problem 3",
      SINGULAR_SYSTEM_BIBTEX,
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 2/(1+power2(x(0)))+sin(x(1)-1)-1;
    f(1) = 2/(1+power2(x(1)))+sin(x(1)-1)-1;
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = -4*x(0)/power2(power2(x(0))+1);
    jac(1) = cos(x(1)-1);
    jac(2) = 0;
    jac(3) = cos(x(1)-1)-4*x(1)/power2(power2(x(1))+1);
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 1;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0.5; x(1) = 0.5; break;
      case 1: x(0) = 2;   x(1) = 2;   break;
      case 2: x(0) = 0.1; x(1) = 0.1; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 3; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP4 : public nonlinearSystem {
public:

  SingularSystemP4()
  : nonlinearSystem(
      "Singular System Problem 4",
      SINGULAR_SYSTEM_BIBTEX,
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 1+tan(2-2*cos(x(0)))-exp(sin(x(0)));
    f(1) = 1+tan(2-2*cos(x(1)))-exp(sin(x(1)));
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = -cos(x(0))*exp(sin(x(0)))+2*sin(x(0))/power2(cos(-2+2*cos(x(0))));
    jac(1) = 0;
    jac(2) = 0;
    jac(3) = -cos(x(1))*exp(sin(x(1)))+2*sin(x(1))/power2(cos(-2+2*cos(x(1))));
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) =    3; x(1) =    0; break;
      case 1: x(0) =    0; x(1) =  0.5; break;
      case 2: x(0) = -0.5; x(1) = -0.5; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 3; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP5 : public nonlinearSystem {
public:

  SingularSystemP5()
  : nonlinearSystem(
      "Singular System Problem 5",
      SINGULAR_SYSTEM_BIBTEX,
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = exp(x(0))+x(1)-1;
    f(1) = exp(x(1))+x(0)-1;
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = exp(x(0));
    jac(1) = 1;
    jac(2) = 1;
    jac(3) = exp(x(1));
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = -0.5; x(1) = -0.5; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP6 : public nonlinearSystem {
public:

  SingularSystemP6()
  : nonlinearSystem(
      "Singular System Problem 6",
      SINGULAR_SYSTEM_BIBTEX,
      3
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    switch ( i ) {
    case 0: return cos(x(0))-9+3*x(0)+8*exp(x(1));
    case 1: return cos(x(1))-9+3*x(1)+8*exp(x(0));
    case 2: return cos(x(2))-x(2)-1;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = cos(x(0))-9+3*x(0)+8*exp(x(1));
    f(1) = cos(x(1))-9+3*x(1)+8*exp(x(0));
    f(2) = cos(x(2))-x(2)-1;
  }

  integer
  jacobianNnz() const override {
    return 5;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(1,0);
    SETIJ(1,1);
    SETIJ(2,2);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = -sin(x(0))+3;
    jac(1) = 8*exp(x(1));
    jac(2) = 8*exp(x(0));
    jac(3) = -sin(x(1))+3;
    jac(4) = -sin(x(2))-1;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
    x(2) = 0;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) =  -1; x(1) =  -1; x(2) =  -1; break;
      case 1: x(0) =   3; x(1) =   3; x(2) =   3; break;
      case 2: x(0) = 0.5; x(1) = 0.5; x(2) = 0.5; break;
      case 3: x(0) =  -3; x(1) =  -3; x(2) =  -3; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 4; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP7 : public nonlinearSystem {
public:

  SingularSystemP7()
  : nonlinearSystem(
      "Singular System Problem 7 (Ishihara, K. 2001)",
      SINGULAR_SYSTEM_BIBTEX,
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = 4*x(0)-2*x(1)+power2(x(0)) - 3;
    f(1) = -2*x(0)+4*x(1)+power2(x(0)) - 3;
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = 4+2*x(0);
    jac(1) = -2;
    jac(2) = -2+2*x(0);
    jac(3) = 4;
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 1;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) =  3; x(1) =    3; break;
      case 1: x(0) =  0; x(1) = -1.5; break;
      case 2: x(0) = -2; x(1) =    3; break;
      case 3: x(0) =  0; x(1) =    2; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 4; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP8 : public nonlinearSystem {
public:

  SingularSystemP8()
  : nonlinearSystem(
      "Singular System Problem 8",
      SINGULAR_SYSTEM_BIBTEX,
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = sqrt(3.0)*power2(x(0)) - power2(x(1));
    f(1) = cos(x(0))-1/(1+power2(x(1)));
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = 2*sqrt(3.0)*x(0);
    jac(1) = -2*x(1);
    jac(2) = -sin(x(0));
    jac(3) = 2*x(1)/power2(power2(x(1))+1);
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) = 0.5; x(1) = 1; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class SingularSystemP9 : public nonlinearSystem {
public:

  SingularSystemP9()
  : nonlinearSystem(
      "Singular System Problem 9",
      SINGULAR_SYSTEM_BIBTEX,
      2
    )
  { }

  real_type
  evalFk( dvec_t const & x, integer i ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(i);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(x(0)) - power2(x(1));
    f(1) = 3*power2(x(0)) - 3*power2(x(1));
  }

  integer
  jacobianNnz() const override
  { return 4; }

  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac(0) = 2*x(0);
    jac(1) = -2*x(1);
    jac(2) = 6*x(0);
    jac(3) = -6*x(1);
  }

  integer
  numExactSolution() const override
  { return 1; }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = 0;
    x(1) = 0;
  }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
      case 0: x(0) =  0.5; x(1) =  0.4; break;
      case 1: x(0) = -0.5; x(1) = -0.4; break;
      case 2: x(0) =  0.3; x(1) = -0.5; break;
      case 3: x(0) =  0.4; x(1) =  0.5; break;
    }
  }

  integer
  numInitialPoint() const override
  { return 4; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
