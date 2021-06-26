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

#define SoniaKrzyworzcka1_BIBTEX \
"@article{Krzyworzcka:1996,\n" \
"  author  = {Sonia Krzyworzcka},\n" \
"  title   = {Extension of the Lanczos and {CGS}\n" \
"             methods to systems of nonlinear equations},\n" \
"  journal = {Journal of Computational and Applied Mathematics},\n" \
"  volume  = {69},\n" \
"  number  = {1},\n" \
"  pages   = {181--190},\n" \
"  year    = {1996},\n" \
"  doi     = {10.1016/0377-0427(95)00032-1}\n" \
"}\n"

class SoniaKrzyworzcka1 : public nonlinearSystem {
public:

  SoniaKrzyworzcka1()
  : nonlinearSystem(
      "Sonia Krzyworzcka example 2",
      SoniaKrzyworzcka1_BIBTEX,
      6
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    switch ( k ) {
    case 0: return -0.75-0.5*x(1)*x(1)*x(3)*x(5) - x(0);
    case 1: return -0.405*exp(1+x(0)*x(1))+1.405 - x(1);
    case 2: return 0.5*x(3)*x(5)-1.5             - x(2);
    case 3: return 0.605*exp(1-x(2)*x(2))+0.395  - x(3);
    case 4: return 0.5*x(1)*x(5)-1.5             - x(4);
    case 5: return x(0)*x(5)                     - x(5);
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override{
    f(0) = -0.75-0.5*x(1)*x(1)*x(3)*x(5) - x(0);
    f(1) = -0.405*exp(1+x(0)*x(1))+1.405 - x(1);
    f(2) = 0.5*x(3)*x(5)-1.5             - x(2);
    f(3) = 0.605*exp(1-x(2)*x(2))+0.395  - x(3);
    f(4) = 0.5*x(1)*x(5)-1.5             - x(4);
    f(5) = x(0)*x(5)                     - x(5);
  }

  integer
  jacobianNnz() const override {
    return 16;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,3);
    SETIJ(0,5);

    SETIJ(1,0);
    SETIJ(1,1);

    SETIJ(2,2);
    SETIJ(2,3);
    SETIJ(2,5);

    SETIJ(3,2);
    SETIJ(3,3);

    SETIJ(4,1);
    SETIJ(4,4);
    SETIJ(4,5);

    SETIJ(5,0);
    SETIJ(5,5);
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    jac(kk++) = -1;
    jac(kk++) = -1.0*x(1)*x(3)*x(5);
    jac(kk++) = -0.5*x(1)*x(1)*x(5);
    jac(kk++) = -0.5*x(1)*x(1)*x(3);
    
    jac(kk++) = -0.405*x(1)*exp(x(0)*x(1)+1);
    jac(kk++) = -0.405*x(0)*exp(x(0)*x(1)+1)-1;

    jac(kk++) = -1;
    jac(kk++) = 0.5*x(5);
    jac(kk++) = 0.5*x(3);
    
    jac(kk++) = -1.210*x(2)*exp(-x(2)*x(2)+1);
    jac(kk++) = -1;

    jac(kk++) = 0.5*x(5);
    jac(kk++) = -1;
    jac(kk++) = 0.5*x(1);
    
    jac(kk++) = x(5);
    jac(kk++) = x(0)-1;
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    x(0) = -1;
    x(1) =  1;
    x(2) = -1;
    x(3) =  1;
    x(4) = -1;
    x(5) =  1;
  }

  integer
  numExactSolution() const override
  { return 1; }
  
  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 0.1;
    x(1) = 0.1;
    x(2) = 0.1;
    x(3) = 0.1;
    x(4) = 0.1;
    x(5) = 0.1;
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

class SoniaKrzyworzcka2 : public nonlinearSystem {
public:

  SoniaKrzyworzcka2()
  : nonlinearSystem(
      "Sonia Krzyworzcka example 3",
      SoniaKrzyworzcka1_BIBTEX,
      10
    )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    for ( integer i = 0; i < n;   ++i ) f(i) = (3-5*x(i))*x(i) + 1;
    for ( integer i = 0; i < n-1; ++i ) f(i) -= 2*x(i+1);
    for ( integer i = 1; i < n;   ++i ) f(i) -= x(i-1);
  }

  integer
  jacobianNnz() const override {
    return 3*n-2;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk
    for ( integer i = 0; i < n;   ++i ) { SETIJ(i,i); }
    for ( integer i = 0; i < n-1; ++i ) { SETIJ(i,i+1); }
    for ( integer i = 1; i < n;   ++i ) { SETIJ(i,i-1); }
    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    integer kk = 0;
    for ( integer i = 0; i < n;   ++i ) jac(kk++) = 3-10*x(i);
    for ( integer i = 0; i < n-1; ++i ) jac(kk++) = -2;
    for ( integer i = 1; i < n;   ++i ) jac(kk++) = -1;
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }
  
  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.fill(-1);
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
