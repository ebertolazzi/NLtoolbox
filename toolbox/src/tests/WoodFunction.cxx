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

class WoodFunction : public nonlinearSystem {
public:

  WoodFunction()
  : nonlinearSystem(
      "Wood function",
      "@book{Colville:1968,\n"
      "  author = {Colville, A.R.},\n"
      "  title  = {A Comparative Study on Nonlinear Programming Codes},\n"
      "  year   = {1968},\n"
      "  notes  = {Rep. 320-2949, New York Scientific Center}\n"
      "}\n\n"
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
      4
    )
  {}
  
  real_type
  t_fun( dvec_t const & x, int_type i ) const {
    switch ( i ) {
    case 0: return sqrt(100.0) * (x(1)-x(0)*x(0));
    case 1: return 1.0 - x(0);
    case 2: return sqrt(90)*(x(3)-x(2)*x(2));
    case 3: return 1.0 - x(2);
    case 4: return sqrt(10)*(x(1)+x(3)-2.0);
    case 5: return sqrt(0.1)* (x(1)-x(3));
    }
    return 0;
  }

  void
  t_grad( dvec_t const & x, int_type i, dvec_t & g ) const {
    g.setZero();
    switch ( i ) {
    case 0:
      g(0) = -20.00*x(0);
      g(1) = 10.00;
    break;
    case 1:
      g(0) = -1;
    break;
    case 2:
      g(2) = -6*sqrt(10)*x(2);
      g(3) = 3*sqrt(10);
    break;
    case 3:
      g(2) = -1;
    break;
    case 4:
      g(1) = g(3) = sqrt(10.0);
    break;
    case 5:
      g(1) = sqrt(0.1);
      g(3) = -sqrt(0.1);
    break;
    }
  }

  void
  t_hess( dvec_t const & x, int_type i, dmat_t & h ) const {
    h.setZero();
    switch ( i ) {
    case 0:
      h(0,0) = -20;
    break;
    case 1:
    break;
    case 2:
      h(2,2) = -6*sqrt(10.0);
    break;
    case 3:
    break;
    case 4:
    break;
    case 5:
    break;
    }
  }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    dvec_t g(4);
    f.setZero();
    for ( int_type i = 0; i < 6; ++i ) {
      real_type t = t_fun( x, i );
      t_grad( x, i, g );
      f += t*g;
    }
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0; // fortran addressing
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    for ( int_type i = 0; i < 6; ++i ) {
      dvec_t g(4);
      dmat_t h(4,4);
      real_type t = t_fun( x, i );
      t_grad( x, i, g );
      t_hess( x, i, h );
      jac[caddr(0,0)] += t*h(0,0)+g(0)*g(0);
      jac[caddr(0,1)] += t*h(0,1)+g(0)*g(1);
      jac[caddr(0,2)] += t*h(0,2)+g(0)*g(2);
      jac[caddr(0,3)] += t*h(0,3)+g(0)*g(3);

      jac[caddr(1,1)] += t*h(1,1)+g(1)*g(1);
      jac[caddr(1,2)] += t*h(1,2)+g(1)*g(2);
      jac[caddr(1,3)] += t*h(1,3)+g(1)*g(3);

      jac[caddr(2,2)] += t*h(2,2)+g(2)*g(2);
      jac[caddr(2,3)] += t*h(2,3)+g(2)*g(3);

      jac[caddr(3,3)] += t*h(3,3)+g(3)*g(3);
    }
    jac[caddr(1,0)] = jac[caddr(0,1)];
    jac[caddr(2,0)] = jac[caddr(0,2)];
    jac[caddr(3,0)] = jac[caddr(0,3)];
    
    jac[caddr(2,1)] = jac[caddr(1,2)];
    jac[caddr(3,1)] = jac[caddr(1,3)];

    jac[caddr(3,2)] = jac[caddr(2,3)];
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 1;
    x(1) = 1;
    x(2) = 1;
    x(3) = 1;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = -3;
    x(1) = -1;
    x(2) = -3;
    x(3) = -1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    // funziona solo con questo limite
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( abs(x(i)) < 100, "x range" );
  }

};
