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

class McKinnon : public nonlinearSystem {
  real_type const tau, theta, phi;
public:

  McKinnon()
  : nonlinearSystem(
      "McKinnon function",
      "@article{McKinnon:1998,\n"
      "  author  = {McKinnon, K.},\n"
      "  title   = {Convergence of the Nelder--Mead Simplex\n"
      "             Method to a Nonstationary Point},\n"
      "  journal = {SIAM Journal on Optimization},\n"
      "  volume  = {9},\n"
      "  number  = {1},\n"
      "  pages   = {148-158},\n"
      "  year    = {1998},\n"
      "  doi     = {10.1137/S1052623496303482},\n"
      "}\n",
      2
    )
  , tau(2.0)
  , theta(6.0)
  , phi(60.0)
  { }

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    if ( x(0) <= 0.0 ) {
      f(0) = theta * tau * phi * pow(-x(0),tau) / x(0);
    } else {
      f(0) = theta * tau * pow(x(0),tau) / x(0);
    }
    f(1) = 1 + (2 * x(1));
  }

  int_type
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
    if ( x(0) <= 0.0 ) {
      real_type t1 = theta * phi;
      real_type t2 = pow(-x(0),tau);
      real_type t3 = tau*tau;
      real_type t5 = x(0)*x(0);
      real_type t6 = 1 / t5;
      jac(0) = t1 * t2 * t3 * t6 - t1 * t2 * tau * t6;
    } else {
      real_type t1 = pow(x(0),tau);
      real_type t2 = theta * t1;
      real_type t3 = tau*tau;
      real_type t4 = x(0)*x(0);
      real_type t5 = 1 / t4;
      jac(0) = t2 * t3 * t5 - t2 * tau * t5;
    }
    jac(1) = 0;
    jac(2) = 0;
    jac(3) = 2;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = -1;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = x(1) = 1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {}

};
