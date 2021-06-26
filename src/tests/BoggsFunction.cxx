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

class BoggsFunction : public nonlinearSystem {
public:
  
  BoggsFunction()
  : nonlinearSystem(
    "Boggs function",
    "@article{Boggs:1971,\n"
    "  author = {Boggs, P.},\n"
    "  title  = {The Solution of Nonlinear Systems of Equations\n"
    "            by A-Stable Integration Techniques},\n"
    "  journal = {SIAM Journal on Numerical Analysis},\n"
    "  volume  = {8},\n"
    "  number  = {4},\n"
    "  pages   = {767--785},\n"
    "  year    = {1971},\n"
    "  doi     = {10.1137/0708071},\n"
    "}\n",
    2
  )
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    switch ( k ) {
     case 0: return power2(x(0))-x(1)+1;
     case 1: return x(0)-cos(m_pi_2*x(1));
    }
    return 0;
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    f(0) = power2(x(0))-x(1)+1;
    f(1) = x(0)-cos(m_pi_2*x(1));
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
    jac(1) = -1;
    jac(2) =  1;
    jac(3) = m_pi_2*sin(m_pi_2*x(1));
  }

  void
  getExactSolution( dvec_t & x, integer  idx ) const override {
    x(0) = 0;
    x(1) = 1;
  }
  
  integer
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 0;
  }

  integer
  numInitialPoint() const override
  { return 1; }


  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
