/*\
 |
 |  Author:
 |    Enrico Bertolazzi
 |    University of Trento
 |    Department of Industrial Engineering
 |    Via Sommarive 9, I-38123, Povo, Trento, Italy
 |    email: enrico.bertolazzi@unitn.it
\*/

class ArtificialTestOfNowakAndWeimann : public nonlinearSystem {
public:
  ArtificialTestOfNowakAndWeimann()
  : nonlinearSystem(
      "Artificial Test of Nowak and Weimann",
      "@techreport{Nowak1991,\n"
      "  author = {U. Nowak and L. Weimann},\n"
      "  title  = {A Family of Newton Codes for Systems of Highly Nonlinear Equations},\n"
      "  number = {Technical Report TR-91-10},\n"
      "  year   = {1991}\n"
      "}\n",
      2
    )
  {}

  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    switch ( k ) {
      case 0: return exp(power2(x(0))+power2(x(1))) - 3;
      case 1: return x(0)+x(1)-sin(3*(x(0)+x(1)));
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = exp(power2(x(0))+power2(x(1))) - 3;
    f(1) = x(0)+x(1)-sin(3*(x(0)+x(1)));
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
  jacobian( dvec_t const & x, dvec_t & vals ) const override {
    real_type tmp1 = 2*exp(power2(x(0))+power2(x(1)));
    vals(0) = x(0)*tmp1;
    vals(1) = x(1)*tmp1;

    real_type tmp2 = 1-3*cos(3*(x(0)+x(1)));
    vals(2) = tmp2;
    vals(3) = tmp2;
  }

  void
  getExactSolution( dvec_t & x, int_type idx ) const override {
  }

  int_type
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = 0.81;
    x(1) = 0.82;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};;
