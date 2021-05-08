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

static
inline
string
ini_msg_MexicanHatFunction( real_type tau ) {
  char msg[1000];
  sprintf( msg, "Mexican Hat Function, tau = %g", tau );
  return string(msg);
}

class MexicanHatFunction : public nonlinearSystem {
  real_type tau;
public:

  MexicanHatFunction( real_type tau_in )
  : nonlinearSystem(
      ini_msg_MexicanHatFunction(tau_in),
      "@article{Grippo:1991,\n"
      "  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n"
      "  title   = {A Class of Nonmonotone Stabilization Methods\n"
      "             in Unconstrained Optimization},\n"
      "  journal = {Numer. Math.},\n"
      "  year    = {1991},\n"
      "  volume  = {59},\n"
      "  number  = {1},\n"
      "  pages   = {779--805},\n"
      "  doi     = {10.1007/BF01385810},\n"
      "}\n",
      2
    )
  , tau(tau_in)
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type t1 = x(0)*x(0);
    real_type t2 = x(1)-t1;
    real_type t3 = t2*t2;
    real_type t6 = power2(1.0-x(0));
    real_type t7 = 10000.0*t3+t6-0.2E-1;
    f(0) = 2*(1-x(0)) + tau*4.0*t7*((1-20000.0*t2)*x(0)-1);
    f(1) = 2*(1-x(1)) + tau*40000.0*t7*t2;
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
    real_type t1  = x(0)*x(0);
    real_type t2  = x(1)-t1;
    real_type t6  = (2-40000*t2)*x(0)-2;
    real_type t7  = t6*t6;
    real_type t8  = t2*t2;
    real_type t11 = power2(1-x(0));
    real_type t12 = 10000*t8+t11-0.02;
    real_type t20 = 20000*(x(1)-t1);
    real_type t25 = 2*t20*t6-80000*t12*x(0);
    real_type t26 = t20*t20;
    jac(0) = -2 + 2*tau*(t7+t12*(120000*t1-40000*x(1)+2));
    jac(1) = tau*t25;
    jac(2) = tau*t25;
    jac(3) = -2 + tau*(2.0*t26+400000000*t8+40000*t11-800);
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type idx ) const override {
    switch ( idx ) {
    case 0:
      x(0) = 0.86;
      x(1) = 0.72;
      break;
    case 1:
      x(0) = 0.85858;
      x(1) = 0.7371534;
      break;
    case 2:
      x(0) = 1.1414204;
      x(1) = 1.3028457;
      break;
    }
  }

  int_type
  numInitialPoint() const override
  { return 3; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
