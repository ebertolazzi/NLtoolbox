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

class PowellQuarticFunction : public nonlinearSystem {

public:

  PowellQuarticFunction()
  : nonlinearSystem(
      "Powell's quartic function",
      "@article{Colville:1970,\n"
      "  author    = {Colville, A. R.},\n"
      "  title     = {A comparative study of nonlinear programming codes},\n"
      "  booktitle = {Proceedings of the {P}rinceton {S}ymposium on\n"
      "               {M}athematical {P}rogramming (1967)},\n"
      "  pages     = {487--501},\n"
      "  publisher = {Princeton Univ. Press, Princeton, N.J.},\n"
      "  year      = {1970},\n"
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
  evalFk( dvec_t const & x_in, int_type k ) const override {
    real_type x = x_in[0];
    real_type y = x_in[1];
    real_type z = x_in[2];
    real_type w = x_in[3];
    
    real_type t4 = x-2.0*y;
    real_type t5 = t4*t4;
    real_type t6 = t5*t4;
    real_type t9 = 10.0*y;
    real_type t10 = 10.0*z;
    real_type t13 = w-z;
    real_type t14 = t13*t13;
    real_type t16 = 40.0*t14*t13;
    switch ( k ) {
      case 0: return 20.0*w+200.0*x+4.0*t6;
      case 1: return t9-t10-8.0*t6;
      case 2: return -t9+t10-t16;
      case 3: return 2.0*w+20.0*x+t16;
    }
    return 0;
  }

  void
  evalF( dvec_t const & x_in, dvec_t & f ) const override {
    real_type x = x_in[0];
    real_type y = x_in[1];
    real_type z = x_in[2];
    real_type w = x_in[3];
    
    //f(0) = 20*w+200*x+4*power3(x-2*y);
    //f(1) = 10*y-10*z-8*power3(x-2*y);
    //f(2) = -10*y+10*z-40*power3(w-z);
    //f(3) = 2*w+20*x+40*power3(w-z);

    real_type t4 = x-2.0*y;
    real_type t5 = t4*t4;
    real_type t6 = t5*t4;
    real_type t9 = 10.0*y;
    real_type t10 = 10.0*z;
    real_type t13 = w-z;
    real_type t14 = t13*t13;
    real_type t16 = 40.0*t14*t13;
    f(0) = 20.0*w+200.0*x+4.0*t6;
    f(1) = t9-t10-8.0*t6;
    f(2) = -t9+t10-t16;
    f(3) = 2.0*w+20.0*x+t16;
  }

  int_type
  jacobianNnz() const override {
    return 12;
  }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    int_type kk = 0;
    #define SETIJ(I,J) ii(kk) = I; jj(kk) = J; ++kk

    SETIJ(0,0);
    SETIJ(0,1);
    SETIJ(0,3);

    SETIJ(1,0);
    SETIJ(1,1);
    SETIJ(1,2);

    SETIJ(2,1);
    SETIJ(2,2);
    SETIJ(2,3);

    SETIJ(3,0);
    SETIJ(3,2);
    SETIJ(3,3);

    #undef SETIJ
  }

  void
  jacobian( dvec_t const & x_in, dvec_t & jac ) const override {
    real_type x = x_in[0];
    real_type y = x_in[1];
    real_type z = x_in[2];
    real_type w = x_in[3];

    real_type t3  = power2(x-2.0*y);
    real_type t6  = 24.0*t3;
    real_type t10 = power2(w-z);
    real_type t11 = 120.0*t10;
    jac(0)  = 200.0+12.0*t3;
    jac(1)  = -t6;
    jac(2)  = 20.0;
    jac(3)  = -t6;
    jac(4)  = 10.0+48.0*t3;
    jac(5)  = -10.0;
    jac(6)  = -10.0;
    jac(7)  = 10.0+t11;
    jac(8)  = -t11;
    jac(9)  = 20.0;
    jac(10) = -t11;
    jac(11) = 2.0+t11;
  }

  void
  getExactSolution( dvec_t & x, int_type ) const override {
    x(0) = 0;
    x(1) = 0;
    x(2) = 0;
    x(3) = 0;
  }

  int_type
  numExactSolution() const override
  { return 1; }

  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) =  3;
    x(1) = -1;
    x(2) =  0;
    x(3) =  1;
  }

  int_type
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override
  {}

};
