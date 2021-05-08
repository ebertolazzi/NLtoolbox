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

class WatsonFunction : public nonlinearSystem {
  real_type t[29][31];
public:

  WatsonFunction()
  : nonlinearSystem(
      "Watson function",
      "@book{brent2013,\n"
      "  author    = {Brent, R.P.},\n"
      "  title     = {Algorithms for Minimization Without Derivatives},\n"
      "  isbn      = {9780486143682},\n"
      "  series    = {Dover Books on Mathematics},\n"
      "  year      = {2013},\n"
      "  publisher = {Dover Publications}\n"
      "}\n\n"
      "@book{kowalik1968methods,\n"
      "  author = {Kowalik, J.S. and Osborne, M.R.},\n"
      "  title  = {Methods for unconstrained optimization problems},\n"
      "  series = {Mathematical Linguistics and Automatic Language Processing},\n"
      "  year   = {1968},\n"
      "  publisher={American Elsevier Pub. Co.}\n"
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
      31
    )
  {
    for ( int_type i = 0; i < 29; ++i ) {
      real_type ti = (i+1)/29.0;
      t[i][0] = 1;
      for ( int_type j = 1; j < 31; ++j ) {
        t[i][j] = t[i][j-1]*ti;
      }
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
    f.setZero();
    for ( int_type i = 0; i < 29; ++i ) {
      real_type const * ti = t[i];
      real_type fi = 0;
      for ( int_type j = 1; j < 31; ++j ) fi += j*x(j)*ti[j-1];
      real_type tmp = 0;
      for ( int_type j = 0; j < 31; ++j ) tmp += x(j)*ti[j];
      f(i) = fi-tmp*tmp-1;
    }
    f(29) = x(0);
    f(30) = x(1) - x(0)*x(0) - 1;
  }

  int_type
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    for ( int_type j = 0; j < n; ++j )
      for ( int_type i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    for ( int_type i = 0; i < 29; ++i ) {
      real_type const * ti = t[i];
      for ( int_type j = 1; j < 31; ++j ) jac(caddr(i,j)) += j*ti[j-1];
      real_type tmp = 0;
      for ( int_type j = 0; j < 31; ++j ) tmp += x(j)*ti[j];
      for ( int_type j = 0; j < 31; ++j ) jac(caddr(i,j)) -= 2*tmp*ti[j];
    }
    jac(caddr(29,0)) = 1;
    jac(caddr(30,0)) = -2*x(0);
    jac(caddr(30,1)) = 1;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    // RES = [ -0.015725, 1.012435, -0.232992, 1.260430, -1.513729, 0.992996 ]';
    // RES = [ -0.000015, 0.999790, 0.014764, 0.146342, 1.000821, -2.617731, 4.104403, -3.143612, 1.052627 ]';
  }

  int_type
  numExactSolution() const override
  { return 0; }
  
  void
  getInitialPoint( dvec_t & x, int_type ) const override
  { x.setZero(); }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    //for (  i = 0; i < n; ++i )
    //  NONLIN_ASSERT( x(i) > -1, "Bad range" );
  }

};
