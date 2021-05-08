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

class EsterificReaction : public nonlinearSystem {

public:

  EsterificReaction()
  : nonlinearSystem(
      "EsterificReaction (example 4)",
      "@book{Luus:2000,\n"
      "  author    = {Luus, Rein},\n"
      "  title     = {Iterative Dynamic Programming},\n"
      "  year      = {2000},\n"
      "  isbn      = {1584881488},\n"
      "  edition   = {1st},\n"
      "  publisher = {CRC Press, Inc.},\n"
      "}\n",
      2
    )
  {}

  virtual
  real_type
  evalFk( dvec_t const & x, int_type k ) const override {
    real_type X  = x(0);
    real_type Y  = x(1);
    real_type t1 = 0.125E1*X;
    real_type t2 = 0.225E1*Y;
    real_type t12 = (X+Y)*(0.6875E1*X+0.7875E1*Y)/(0.5797E1*X+0.6797E1*Y);
    real_type t16 = pow(t12-Y-X,2.0);
    real_type t20 = pow(0.11526E1*t12-t1-t2,2.0);
    switch ( k ) {
      case 0: return (t1+t2-0.1054E1*t12)*t16-0.20564E4*t20;
      case 1: return 0.55E1*(-0.1433314531E3*t12+0.2801659319E3*Y+0.7647126E2+0.1556548821E3*X)*
                            (-0.1043629412E3*t12+0.2031305325E3*Y+0.1134099326E3*X+0.61177E2)
                             - (0.1422774531E3*t12-0.2789159319E3*Y-0.1544048821E3*X-0.7585949E2)
                             * (0.11526E3*t12-0.225E3*Y-0.125E3*X-0.61177E2);
    }
    return 0;
  }

  virtual
  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type X = x(0);
    real_type Y = x(1);
    real_type t1 = 0.125E1*X;
    real_type t2 = 0.225E1*Y;
    real_type t12 = (X+Y)*(0.6875E1*X+0.7875E1*Y)/(0.5797E1*X+0.6797E1*Y);
    real_type t16 = pow(t12-Y-X,2.0);
    real_type t20 = pow(0.11526E1*t12-t1-t2,2.0);
    f(0) = (t1+t2-0.1054E1*t12)*t16-0.20564E4*t20;
    f(1) = 0.55E1*(-0.1433314531E3*t12+0.2801659319E3*Y+0.7647126E2+0.1556548821E3*X)*
                  (-0.1043629412E3*t12+0.2031305325E3*Y+0.1134099326E3*X+0.61177E2)
                   - (0.1422774531E3*t12-0.2789159319E3*Y-0.1544048821E3*X-0.7585949E2)
                   * (0.11526E3*t12-0.225E3*Y-0.125E3*X-0.61177E2);
  }

  virtual
  int_type
  jacobianNnz() const override
  { return 4; }

  virtual
  void
  jacobianPattern( ivec_t & i, ivec_t & j ) const override {
    i(0) = 0; j(0) = 0;
    i(1) = 0; j(1) = 1;
    i(2) = 1; j(2) = 0;
    i(3) = 1; j(3) = 1;
  }

  virtual
  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type X = x(0);
    real_type Y = x(1);
    real_type t3 = 0.6875E1*X+0.7875E1*Y;
    real_type t6 = 0.5797E1*X+0.6797E1*Y;
    real_type t7 = 1/t6;
    real_type t8 = t3*t7;
    real_type t9 = 0.1054E1*t8;
    real_type t10 = X+Y;
    real_type t11 = t10*t7;
    real_type t13 = t10*t3;
    real_type t14 = t6*t6;
    real_type t16 = t13/t14;
    real_type t19 = t13*t7;
    real_type t20 = t19-Y-X;
    real_type t21 = t20*t20;
    real_type t23 = 0.125E1*X;
    real_type t24 = 0.225E1*Y;
    real_type t27 = (t23+t24-0.1054E1*t19)*t20;
    real_type t34 = 0.11526E1*t19-t23-t24;
    real_type t35 = 0.11526E1*t8;
    real_type t57 = 0.1433314531E3*t8;
    real_type t64 = -0.1043629412E3*t19+0.2031305325E3*Y+0.1134099326E3*X+0.61177E2;
    real_type t70 = -0.1433314531E3*t19+0.2801659319E3*Y+0.7647126E2+0.1556548821E3*X;
    real_type t71 = 0.1043629412E3*t8;
    real_type t77 = 0.1422774531E3*t8;
    real_type t84 = 0.11526E3*t19-0.225E3*Y-0.125E3*X-0.61177E2;
    real_type t89 = 0.1422774531E3*t19-0.2789159319E3*Y-0.1544048821E3*X-0.7585949E2;
    real_type t90 = 0.11526E3*t8;
    jac(0) = (0.125E1-t9-0.724625E1*t11+0.6110038E1*t16)*t21
                    +2.0*t27*(t8+0.6875E1*t11-0.5797E1*t16-1.0)
                    -0.41128E4*t34*(t35+0.7924125E1*t11-0.66816222E1*t16-0.125E1);
    jac(1) = (0.225E1-t9-0.830025E1*t11+0.7164038E1*t16)*t21
                    +2.0*t27*(t8+0.7875E1*t11-0.6797E1*t16-1.0)
                    -0.41128E4*t34*(t35+0.9076725E1*t11-0.78342222E1*t16-0.225E1);
    jac(2) = 0.55E1*(-t57-0.9854037401E3*t11+0.8308924336E3*t16+0.1556548821E3)*t64
                    + 0.55E1*t70*(-t71-0.7174952208E3*t11+0.6049919701E3*t16+0.1134099326E3)
                    - (t77+0.9781574901E3*t11-0.8247823956E3*t16-0.1544048821E3)*t84
                    -t89*(t90+0.7924125E3*t11-0.66816222E3*t16-0.125E3);
    jac(3) = 0.55E1*(-t57-0.1128735193E4*t11+0.9742238867E3*t16+0.2801659319E3)*t64
                    + 0.55E1*t70*(-t71-0.821858162E3*t11+0.7093549113E3*t16+0.2031305325E3)
                    - (t77+0.1120434943E4*t11-0.9670598487E3*t16-0.2789159319E3)*t84
                    - t89*(t90+0.9076725E3*t11-0.78342222E3*t16-0.225E3);
  }

  virtual
  int_type
  numExactSolution() const override
  { return 4; }

  virtual
  void
  getExactSolution( dvec_t & x ,int_type idx ) const override {
    static real_type const s0[2] = {61.15373203, 6.976291975};
    static real_type const s1[2] = {67.09266777, 7.614805557};
    static real_type const s2[2] = {20.85605925, -24.33016974};
    static real_type const s3[2] = {-33.28012847, 28.52767665};
    switch ( idx ) {
      case 0: x(0) = s0[0]; x(1) = s0[1]; break;
      case 1: x(0) = s1[0]; x(1) = s1[1]; break;
      case 2: x(0) = s2[0]; x(1) = s2[1]; break;
      case 3: x(0) = s3[0]; x(1) = s3[1]; break;
    }
  }

  virtual
  void
  getInitialPoint( dvec_t & x, int_type ) const override {
    x(0) = x(1) = 20;
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
