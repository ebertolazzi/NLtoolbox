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

#define HAS_BIBTEX \
"@article{Grippo:1991,\n" \
"  author  = {Grippo, L. and Lampariello, F. and Lucidi, S.},\n" \
"  title   = {A Class of Nonmonotone Stabilization Methods\n" \
"             in Unconstrained Optimization},\n" \
"  journal = {Numer. Math.},\n" \
"  year    = {1991},\n" \
"  volume  = {59},\n" \
"  number  = {1},\n" \
"  pages   = {779--805},\n" \
"  doi     = {10.1007/BF01385810},\n" \
"}\n"

class HAS64 : public nonlinearSystem {
  typedef pair<integer,integer> INDEX;
  mutable map<INDEX,real_type> jac_idx_vals;
  real_type tau;
public:

  HAS64( real_type tau_in)
  : nonlinearSystem(
      fmt::format( "HAS 64, tau = {}", tau_in ),
      HAS_BIBTEX,
      7
    )
  , tau(tau_in)
  {
    jac_idx_vals.clear();
    jac_idx_vals[INDEX(0,0)] = 1;
    jac_idx_vals[INDEX(0,4)] = 1;
    jac_idx_vals[INDEX(1,1)] = 1;
    jac_idx_vals[INDEX(1,5)] = 1;
    jac_idx_vals[INDEX(2,2)] = 1;
    jac_idx_vals[INDEX(2,6)] = 1;
    jac_idx_vals[INDEX(4,4)] = 1;
    jac_idx_vals[INDEX(4,0)] = 1;
    jac_idx_vals[INDEX(5,5)] = 1;
    jac_idx_vals[INDEX(5,1)] = 1;
    jac_idx_vals[INDEX(6,6)] = 1;
    jac_idx_vals[INDEX(6,2)] = 1;

    jac_idx_vals[INDEX(0,0)] = 1;
    jac_idx_vals[INDEX(0,1)] = 1;
    jac_idx_vals[INDEX(0,2)] = 1;
    jac_idx_vals[INDEX(0,3)] = 1;

    jac_idx_vals[INDEX(1,0)] = 1;
    jac_idx_vals[INDEX(1,1)] = 1;
    jac_idx_vals[INDEX(1,2)] = 1;
    jac_idx_vals[INDEX(1,3)] = 1;

    jac_idx_vals[INDEX(2,0)] = 1;
    jac_idx_vals[INDEX(2,1)] = 1;
    jac_idx_vals[INDEX(2,2)] = 1;
    jac_idx_vals[INDEX(2,3)] = 1;

    jac_idx_vals[INDEX(3,0)] = 1;
    jac_idx_vals[INDEX(3,1)] = 1;
    jac_idx_vals[INDEX(3,2)] = 1;
    jac_idx_vals[INDEX(3,3)] = 1;

    jac_idx_vals[INDEX(0,0)] = 1;
    jac_idx_vals[INDEX(1,1)] = 1;
    jac_idx_vals[INDEX(2,2)] = 1;

  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type x0 = x(0);
    real_type x1 = x(1);
    real_type x2 = x(2);
    real_type x3 = x(3);
    real_type x4 = x(4);
    real_type x5 = x(5);
    real_type x6 = x(6);
    if ( x0 <= 0 || x1 <= 0 || x2 <= 0 ) {
      f(0) = f(1) = f(2) = f(3) = f(4) = f(5) = f(6) = nan("HAS64");
      return;
    }

    real_type x0x0 = x0*x0;
    real_type x1x1 = x1*x1;
    real_type x2x2 = x2*x2;
    real_type x3x3 = x3*x3;
    real_type x4x4 = x4*x4;
    real_type x5x5 = x5*x5;
    real_type x6x6 = x6*x6;

    f(0) = 2*(x0-x4x4)-2E-5;
    f(1) = 2*(x1-x5x5)-2E-5;
    f(2) = 2*(x2-x6x6)-2E-5;
    f(3) = 0.0;
    f(4) = 4*(x4x4-x0+1E-5)*x4;
    f(5) = 4*(x5x5-x1+1E-5)*x5;
    f(6) = 4*(x6x6-x2+1E-5)*x6;

    real_type tmp = 1 - 4/x0 - 32/x1 - 120/x2 - x3x3;

    f(0) += tmp*(8/x0x0);
    f(1) += tmp*(64/x1x1);
    f(2) += tmp*(240/x2x2);
    f(3) -= 4*tmp*x3;

    f *= tau;

    f(0) += 5;
    f(1) += 20;
    f(2) += 10;

    f(0) -= 50000/x0x0;
    f(1) -= 72000/x1x1;
    f(2) -= 144000/x2x2;
  }

  integer
  jacobianNnz() const override
  { return integer(jac_idx_vals.size()); }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    map<pair<integer,integer>,real_type>::const_iterator it = jac_idx_vals.begin();
    integer kk = 0;
    for (; it != jac_idx_vals.end(); ++it ) {
      ii(kk) = it->first.first;
      jj(kk) = it->first.second;
      ++kk;
    }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    map<pair<integer,integer>,real_type>::iterator it = jac_idx_vals.begin();
    for (; it != jac_idx_vals.end(); ++it ) it->second = 0;

    real_type x0 = x(0);
    real_type x1 = x(1);
    real_type x2 = x(2);
    real_type x3 = x(3);
    real_type x4 = x(4);
    real_type x5 = x(5);
    real_type x6 = x(6);

    real_type x0x0 = x0*x0;
    real_type x1x1 = x1*x1;
    real_type x2x2 = x2*x2;
    real_type x3x3 = x3*x3;
    real_type x4x4 = x4*x4;
    real_type x5x5 = x5*x5;
    real_type x6x6 = x6*x6;

    //f(0) += 5-50000/power2(x(0));
    //f(1) += 20-72000/power2(x(1));
    //f(2) += 10-144000/power2(x(2));
    jac_idx_vals[INDEX(0,0)] += 100000/(x0*x0x0);
    jac_idx_vals[INDEX(1,1)] += 144000/(x1*x1x1);
    jac_idx_vals[INDEX(2,2)] += 288000/(x2*x2x2);

    jac_idx_vals[INDEX(0,0)] += 2*tau;
    jac_idx_vals[INDEX(0,4)] -= 4*x4*tau;

    jac_idx_vals[INDEX(1,1)] += 2*tau;
    jac_idx_vals[INDEX(1,5)] -= 4*x5*tau;

    jac_idx_vals[INDEX(2,2)] += 2*tau;
    jac_idx_vals[INDEX(2,6)] -= 4*x6*tau;

    jac_idx_vals[INDEX(4,0)] -= 4*x4*tau;
    jac_idx_vals[INDEX(4,4)] += 4*(3*x4x4-x0+0.00001)*tau;

    jac_idx_vals[INDEX(5,1)] -= 4*x5*tau;
    jac_idx_vals[INDEX(5,5)] += 4*(3*x5x5-x1+0.00001)*tau;

    jac_idx_vals[INDEX(6,2)] -= 4*x6*tau;
    jac_idx_vals[INDEX(6,6)] += 4*(3*x6x6-x2+0.00001)*tau;

    real_type tmp   = 1 - 4/x0 - 32/x1 -120/x2 - x3x3;
    real_type tmp_0 = 4/x0x0;
    real_type tmp_1 = 32/x1x1;
    real_type tmp_2 = 120/x2x2;
    real_type tmp_3 = -2*x3;

    real_type tt1 = 8*tau/x0x0;
    jac_idx_vals[INDEX(0,0)] += tt1*(tmp_0-2*tmp/x0);
    jac_idx_vals[INDEX(0,1)] += tt1*tmp_1;
    jac_idx_vals[INDEX(0,2)] += tt1*tmp_2;
    jac_idx_vals[INDEX(0,3)] += tt1*tmp_3;

    real_type tt2 = 64*tau/x1x1;
    jac_idx_vals[INDEX(1,0)] += tt2*tmp_0;
    jac_idx_vals[INDEX(1,1)] += tt2*(tmp_1-2*tmp/x1);
    jac_idx_vals[INDEX(1,2)] += tt2*tmp_2;
    jac_idx_vals[INDEX(1,3)] += tt2*tmp_3;

    real_type tt3 = 240*tau/x2x2;
    jac_idx_vals[INDEX(2,0)] += tt3*tmp_0;
    jac_idx_vals[INDEX(2,1)] += tt3*tmp_1;
    jac_idx_vals[INDEX(2,2)] += tt3*(tmp_2-2*tmp/x2);
    jac_idx_vals[INDEX(2,3)] += tt3*tmp_3;

    real_type tt4 = 4*tau;
    jac_idx_vals[INDEX(3,0)] -= tt4*(tmp_0*x3);
    jac_idx_vals[INDEX(3,1)] -= tt4*(tmp_1*x3);
    jac_idx_vals[INDEX(3,2)] -= tt4*(tmp_2*x3);
    jac_idx_vals[INDEX(3,3)] -= tt4*(tmp-2*x3x3);

    integer kk = 0;
    for ( it = jac_idx_vals.begin(); it != jac_idx_vals.end(); ++it )
      jac(kk++) = it->second;
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
    if ( tau == 1e2 ) {
      x << 6.6818141799970216926, 9.0745104956268746744, 12.995355396076109065, 0, 0, 0, 0;
    } if ( tau == 1e4 ) {
      x << 3.9001579154194764056, 7.7092485438427308755, 11.962907630149328334, 0, 0, 0, 0;
    } if ( tau == 1e6 ) {
      x << 3.847300847063338, 7.694589110936490, 11.954435422989896, 0, 0, 0, 0;
    } if ( tau == 1e8 ) {
      x << 3.8473008470633382174, 7.6945891109364897632, 11.954435422989895821, 0, 0, 0, 0;
    } if ( tau == 1e10 ) {
      x << 3.8472955440686861651, 7.6945876623010364843, 11.954434592374304256, 0, 0, 0, 0;
    }
  }

  integer
  numExactSolution() const override {
    if ( tau == 1e2 || tau == 1e4 || tau == 1e6 || tau == 1e8 || tau == 1e10 ) return 1;
    return 0;
  }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x(0) = 1;
    x(1) = 1;
    x(2) = 1;
    x(3) = -10;
    x(4) = -10;
    x(5) = -10;
    x(6) = -10;
  }

  integer
  numInitialPoint() const override
  { return 1; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    UTILS_ASSERT( x(0) > 0, "checkIfAdmissible x(0) = {} must be > 0", x(0) );
    UTILS_ASSERT( x(1) > 0, "checkIfAdmissible x(1) = {} must be > 0", x(1) );
    UTILS_ASSERT( x(2) > 0, "checkIfAdmissible x(2) = {} must be > 0", x(2) );
    // UTILS_ASSERT( x(3) < 0, "Bad range" );
    // UTILS_ASSERT( x(4) < 0, "Bad range" );
    // UTILS_ASSERT( x(5) < 0, "Bad range" );
    // UTILS_ASSERT( x(6) < 0, "Bad range" );
  }

  void
  boundingBox( dvec_t & L, dvec_t & U ) const override {
    U.fill(real_max);
    L.fill(-real_max);
    L(0) = L(1) = L(2) = 0;
  }

};

/*\
 | - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
\*/

class HAS93 : public nonlinearSystem {
  real_type tau;
public:

  HAS93( real_type tau_in)
  : nonlinearSystem(
      fmt::format( "HAS 93, tau = {}", tau_in),
      HAS_BIBTEX,
      14
    )
  , tau(tau_in)
  {}

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
   {
      real_type t4 = x(3)*x(4)*x(5);
      real_type t7 = x(6)*x(6);
      real_type t8 = 0.1E-2*x(0)*x(1)*x(2)*t4-0.207E1-t7;
      real_type t13 = x(0)*x(3);
      real_type t14 = x(4)*x(4);
      real_type t15 = x(0)+x(1)+x(2);
      real_type t19 = x(1)*x(2);
      real_type t20 = x(5)*x(5);
      real_type t22 = x(0)+0.157E1*x(1)+x(3);
      real_type t26 = x(7)*x(7);
      real_type t27 = 1.0-0.62E-3*t13*t14*t15-0.58E-3*t19*t20*t22-t26;
      real_type t32 = 0.62E-3*t13*t14;
      real_type t33 = t19*t20;
      real_type t34 = 0.58E-3*t33;
      real_type t41 = t8*x(0);
      real_type t55 = t41*x(1);
      real_type t80 = x(2)*x(3);
      f(0)  = 0.2E-2*t8*x(1)*x(2)*t4+2.0*t27*(-0.62E-3*x(3)*t14*t15-t32-t34)+2.0*x(0)-2.0*x(8);
      f(1)  = 0.2E-2*t41*x(2)*t4+2.0*t27*(-t32-0.58E-3*x(2)*t20*t22-0.9106E-3*t33)+2.0*x(1)-2.0*x(9);
      f(2)  = 0.2E-2*t55*t4+2.0*t27*(-t32-0.58E-3*x(1)*t20*t22)+2.0*x(2)-2.0*x(10);
      f(3)  = 0.2E-2*t55*x(2)*x(4)*x(5)+2.0*t27*(-0.62E-3*x(0)*t14*t15-t34)+2.0*x(3)-2.0*x(11);
      f(4)  = 0.2E-2*t55*t80*x(5)-0.248E-2*t27*x(0)*x(3)*t15*x(4)+2.0*x(4)-2.0*x(12);
      f(5)  = 0.2E-2*t55*t80*x(4)-0.232E-2*t27*x(1)*x(2)*t22*x(5)+2.0*x(5)-2.0*x(13);
      f(6)  = -4.0*t8*x(6);
      f(7)  = -4.0*t27*x(7);
      f(8)  = -2.0*x(0)+2.0*x(8);
      f(9)  = -2.0*x(1)+2.0*x(9);
      f(10) = -2.0*x(2)+2.0*x(10);
      f(11) = -2.0*x(3)+2.0*x(11);
      f(12) = -2.0*x(4)+2.0*x(12);
      f(13) = -2.0*x(5)+2.0*x(13);
    }

    for ( integer i = 0; i < n; ++i ) f(i) *= tau;

    real_type t2  = x(4)*x(4);
    real_type t13 = x(2)*x(1);
    real_type t14 = x(5)*x(5);
    real_type t16 = 0.187E-1+0.437E-1*t14;
    real_type t17 = t13*t16;
    real_type t30 = x(3)*x(0);
    real_type t32 = 0.204E-1+0.607E-1*t2;
    real_type t33 = t30*t32;
    real_type t36 = x(0)+0.157E1*x(1)+x(3);
    real_type t40 = x(0)+x(1)+x(2);
    f(0) += (0.408E-1*x(0)+0.1214*x(0)*t2+0.204E-1*x(1)+0.607E-1*x(1)*t2+0.204E-1*x(2)+0.607E-1*x(2)*t2)*x(3)+t17;
    f(1) += (0.187E-1*x(0)+0.437E-1*x(0)*t14+0.58718E-1*x(1)+0.137218*x(1)*t14+0.187E-1*x(3)+0.437E-1*x(3)*t14)*x(2)+t33;
    f(2) += x(1)*t36*t16+t33;
    f(3) += x(0)*t40*t32+t17;
    f(4) += 0.1214*t30*t40*x(4);
    f(5) += 0.874E-1*t13*t36*x(5);
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0; // fortran address
    for ( integer j = 0; j < n; ++j )
      for ( integer i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    jac.setZero();
    {
      real_type t1 = x(1)*x(1);
      real_type t2 = x(2)*x(2);
      real_type t3 = t1*t2;
      real_type t4 = x(3)*x(3);
      real_type t5 = x(4)*x(4);
      real_type t6 = t4*t5;
      real_type t7 = x(5)*x(5);
      real_type t8 = t6*t7;
      real_type t11 = x(3)*t5;
      real_type t12 = x(0)+x(1)+x(2);
      real_type t15 = x(0)*x(3);
      real_type t17 = 0.62E-3*t15*t5;
      real_type t18 = x(1)*x(2);
      real_type t19 = t18*t7;
      real_type t20 = 0.58E-3*t19;
      real_type t21 = -0.62E-3*t11*t12-t17-t20;
      real_type t22 = t21*t21;
      real_type t24 = t5*t5;
      real_type t27 = x(2)*t7;
      real_type t28 = t27*t5;
      real_type t32 = x(0)*x(0);
      real_type t40 = t7*x(0);
      real_type t46 = x(7)*x(7);
      real_type t54 = t5*t7;
      real_type t60 = x(3)*x(4);
      real_type t61 = t60*x(5);
      real_type t62 = x(0)*x(1)*x(2)*t61;
      real_type t64 = x(6)*x(6);
      real_type t65 = 0.1E-2*t62-0.207E1-t64;
      real_type t70 = x(0)+0.157E1*x(1)+x(3);
      real_type t74 = -t17-0.58E-3*t27*t70-0.9106E-3*t19;
      real_type t77 = t5*t12;
      real_type t78 = t15*t77;
      real_type t80 = t7*t70;
      real_type t81 = t18*t80;
      real_type t83 = 1.0-0.62E-3*t78-0.58E-3*t81-t46;
      real_type t84 = 0.62E-3*t11;
      real_type t85 = 0.58E-3*t27;
      real_type t89 = 0.2E-5*x(0)*t2*t4*t54*x(1)+0.2E-2*t65*x(2)*t61+2.0*t74*t21+2.0*t83*(-t84-t85);
      real_type t90 = x(0)*t1;
      real_type t94 = t65*x(1);
      real_type t97 = t7*x(1);
      real_type t100 = -t17-0.58E-3*t97*t70;
      real_type t103 = 0.58E-3*t97;
      real_type t107 = 0.2E-5*t90*t4*t28+0.2E-2*t94*t61+2.0*t100*t21+2.0*t83*(-t84-t103);
      real_type t108 = t90*t2;
      real_type t109 = t54*x(3);
      real_type t113 = x(2)*x(4)*x(5);
      real_type t116 = x(0)*t5;
      real_type t119 = -0.62E-3*t116*t12-t20;
      real_type t123 = 0.62E-3*t116;
      real_type t127 = 0.2E-5*t108*t109+0.2E-2*t94*t113+2.0*t119*t21+2.0*t83*(-0.62E-3*t77-t123);
      real_type t129 = t4*t7*x(4);
      real_type t131 = 0.2E-5*t108*t129;
      real_type t132 = x(2)*x(3);
      real_type t133 = t132*x(5);
      real_type t135 = 0.2E-2*t94*t133;
      real_type t136 = t12*x(4);
      real_type t139 = 0.248E-2*t15*t136*t21;
      real_type t140 = t60*t12;
      real_type t148 = t6*x(5);
      real_type t151 = t132*x(4);
      real_type t154 = t70*x(5);
      real_type t158 = t83*x(1);
      real_type t159 = x(2)*x(5);
      real_type t160 = t158*t159;
      real_type t161 = 0.232E-2*t160;
      real_type t162 = 0.2E-5*t108*t148+0.2E-2*t94*t151-0.232E-2*t18*t154*t21-t161;
      real_type t166 = 0.4E-2*x(6)*x(1)*x(2)*t61;
      real_type t173 = x(7)*x(1);
      real_type t176 = x(7)*(0.496E-2*x(0)+0.248E-2*x(1)+0.248E-2*x(2))*t11+0.232E-2*t173*t27;
      real_type t180 = t74*t74;
      real_type t182 = t7*t7;
      real_type t206 = t32*x(1);
      real_type t210 = t65*x(0);
      real_type t220 = 0.2E-5*t206*t4*t28+0.2E-2*t210*t61+2.0*t100*t74+2.0*t83*(-0.58E-3*t80-0.9106E-3*t97);
      real_type t221 = t206*t2;
      real_type t231 = 0.2E-5*t221*t109+0.2E-2*t210*t113+2.0*t119*t74+2.0*t83*(-t123-t85);
      real_type t239 = t83*x(0);
      real_type t241 = 0.248E-2*t239*t60;
      real_type t242 = 0.2E-5*t221*t129+0.2E-2*t210*t133-0.248E-2*t15*t136*t74-t241;
      real_type t244 = 0.2E-5*t221*t148;
      real_type t246 = 0.2E-2*t210*t151;
      real_type t249 = 0.232E-2*t18*t154*t74;
      real_type t250 = t159*t70;
      real_type t252 = t18*x(5);
      real_type t258 = x(6)*x(0);
      real_type t261 = 0.4E-2*t258*x(2)*t61;
      real_type t272 = 0.232E-2*x(7)*(0.1E1*x(0)+0.314E1*x(1)+0.1E1*x(3))*t27+0.248E-2*x(7)*x(3)*t116;
      real_type t273 = t32*t1;
      real_type t276 = t100*t100;
      real_type t279 = t273*x(2);
      real_type t291 = 0.2E-5*t279*t109+0.2E-2*t210*x(1)*x(4)*x(5)+2.0*t119*t100+2.0*t83*(-t123-t103);
      real_type t294 = x(1)*x(3);
      real_type t301 = 0.2E-5*t279*t129+0.2E-2*t210*t294*x(5)-0.248E-2*t15*t136*t100-t241;
      real_type t312 = 0.2E-5*t279*t148+0.2E-2*t210*t294*x(4)-0.232E-2*t18*t154*t100-0.232E-2*t158*t154;
      real_type t313 = t258*x(1);
      real_type t315 = 0.4E-2*t313*t61;
      real_type t317 = 4.0*x(7)*t100;
      real_type t322 = t119*t119;
      real_type t325 = t273*t2;
      real_type t337 = 0.2E-5*t325*x(3)*t7*x(4)+0.2E-2*t210*t252-0.248E-2*t15*t136*t119-0.248E-2*t239*t136;
      real_type t347 = 0.2E-5*t325*t11*x(5)+0.2E-2*t210*t18*x(4)-0.232E-2*t18*t154*t119-t161;
      real_type t349 = 0.4E-2*t313*t113;
      real_type t351 = 4.0*x(7)*t119;
      real_type t356 = t2*t4;
      real_type t361 = t12*t12;
      real_type t380 = 0.2E-5*t325*t4*x(4)*x(5)+0.2E-2*t210*t18*x(3)+0.28768E-5*t18*t154*t15*t136;
      real_type t382 = 0.4E-2*t313*t133;
      real_type t385 = 0.496E-2*x(7)*x(0)*t140;
      real_type t394 = t70*t70;
      real_type t403 = 0.4E-2*t313*t151;
      real_type t405 = 0.464E-2*t173*t250;

      jac[caddr(0,0)]   = 0.2E-5*t3*t8+2.0*t22+((0.15376E-5*t24*x(0)+0.14384E-5*t28)*x(1)+0.15376E-5*t24*t32+0.15376E-5*t24*x(2)*x(0))*t4+0.248E-2*(-1.0+0.58E-3*t18*t40+0.9106E-3*t1*x(2)*t7+t46)*t5*x(3)+2.0;
      jac[caddr(0,1)]   = t89;
      jac[caddr(0,2)]   = t107;
      jac[caddr(0,3)]   = t127;
      jac[caddr(0,4)]   = t131+t135-t139+2.0*t83*(-0.124E-2*t140-0.124E-2*t15*x(4));
      jac[caddr(0,5)]   = t162;
      jac[caddr(0,6)]   = -t166;
      jac[caddr(0,7)]   = t176;
      jac[caddr(0,8)]   = -2.0;

      jac[caddr(1,0)]   = t89;
      jac[caddr(1,1)]   = 0.2E-5*t32*t2*t8+2.0*t180+(0.331676944E-5*t1*t182+(0.2112592E-5*t182*x(0)+0.2112592E-5*t182*x(3))*x(1)+0.2258288E-5*t40*t11)*t2+0.36424E-2*(-1.0+0.62E-3*x(3)*t32*t5+0.62E-3*t15*t5*x(1)+t46)*t7*x(2)+2.0;
      jac[caddr(1,2)]   = t220;
      jac[caddr(1,3)]   = t231;
      jac[caddr(1,4)]   = t242;
      jac[caddr(1,5)]   = t244+t246-t249+2.0*t83*(-0.116E-2*t250-0.18212E-2*t252);
      jac[caddr(1,6)]   = -t261;
      jac[caddr(1,7)]   = t272;
      jac[caddr(1,9)]   = -2.0;

      jac[caddr(2,0)]   = t107;
      jac[caddr(2,1)]   = t220;
      jac[caddr(2,2)]   = 0.2E-5*t273*t8+2.0*t276+2.0;
      jac[caddr(2,3)]   = t291;
      jac[caddr(2,4)]   = t301;
      jac[caddr(2,5)]   = t312;
      jac[caddr(2,6)]   = -t315;
      jac[caddr(2,7)]   = -t317;
      jac[caddr(2,10)]  = -2.0;

      jac[caddr(3,0)]   = t127;
      jac[caddr(3,1)]   = t231;
      jac[caddr(3,2)]   = t291;
      jac[caddr(3,3)]   = 0.2E-5*t273*t2*t5*t7+2.0*t322+2.0;
      jac[caddr(3,4)]   = t337;
      jac[caddr(3,5)]   = t347;
      jac[caddr(3,6)]   = -t349;
      jac[caddr(3,7)]   = -t351;
      jac[caddr(3,11)]  = -2.0;

      jac[caddr(4,0)]   = t131+t135-t139-0.248E-2*t83*x(3)*t136-t241;
      jac[caddr(4,1)]   = t242;
      jac[caddr(4,2)]   = t301;
      jac[caddr(4,3)]   = t337;
      jac[caddr(4,4)]   = 0.2E-5*t273*t356*t7+0.30752E-5*t32*t4*t361*t5-0.248E-2*t239*x(3)*t12+2.0;
      jac[caddr(4,5)]   = t380;
      jac[caddr(4,6)]   = -t382;
      jac[caddr(4,7)]   = t385;
      jac[caddr(4,12)]  = -2.0;

      jac[caddr(5,0)]   = t162;
      jac[caddr(5,1)]   = t244+t246-t249-0.232E-2*t83*x(2)*t154-0.36424E-2*t160;
      jac[caddr(5,2)]   = t312;
      jac[caddr(5,3)]   = t347;
      jac[caddr(5,4)]   = t380;
      jac[caddr(5,5)]   = 0.2E-5*t273*t356*t5+0.26912E-5*t3*t394*t7-0.232E-2*t158*x(2)*t70+2.0;
      jac[caddr(5,6)]   = -t403;
      jac[caddr(5,7)]   = t405;
      jac[caddr(5,13)]  = -2.0;

      jac[caddr(6,0)]   = -t166;
      jac[caddr(6,1)]   = -t261;
      jac[caddr(6,2)]   = -t315;
      jac[caddr(6,3)]   = -t349;
      jac[caddr(6,4)]   = -t382;
      jac[caddr(6,5)]   = -t403;
      jac[caddr(6,6)]   = 12.0*t64-0.4E-2*t62+0.828E1;

      jac[caddr(7,0)]   = t176;
      jac[caddr(7,1)]   = t272;
      jac[caddr(7,2)]   = -t317;
      jac[caddr(7,3)]   = -t351;
      jac[caddr(7,4)]   = t385;
      jac[caddr(7,5)]   = t405;
      jac[caddr(7,7)]   = 12.0*t46-4.0+0.248E-2*t78+0.232E-2*t81;

      jac[caddr(8,0)]   = -2.0;
      jac[caddr(8,8)]   = 2.0;

      jac[caddr(9,1)]   = -2.0;
      jac[caddr(9,9)]   = 2.0;

      jac[caddr(10,2)]  = -2.0;
      jac[caddr(10,10)] = 2.0;

      jac[caddr(11,3)]  = -2.0;
      jac[caddr(11,11)] = 2.0;

      jac[caddr(12,4)]  = -2.0;
      jac[caddr(12,12)] = 2.0;

      jac[caddr(13,5)]  = -2.0;
      jac[caddr(13,13)] = 2.0;
    }

    jac *= tau;

    {
      real_type t1 = x(4)*x(4);
      real_type t2 = x(3)*t1;
      real_type t6 = 0.204E-1*x(3);
      real_type t7 = 0.607E-1*t2;
      real_type t8 = 0.187E-1*x(2);
      real_type t9 = x(5)*x(5);
      real_type t10 = x(2)*t9;
      real_type t11 = 0.437E-1*t10;
      real_type t12 = t6+t7+t8+t11;
      real_type t13 = 0.187E-1*x(1);
      real_type t15 = 0.437E-1*x(1)*t9;
      real_type t16 = t6+t7+t13+t15;
      real_type t31 = x(4)*(0.2428*x(0)+0.1214*x(1)+0.1214*x(2))*x(3);
      real_type t32 = x(1)*x(2);
      real_type t34 = 0.874E-1*t32*x(5);
      real_type t47 = 0.204E-1*x(0);
      real_type t49 = 0.607E-1*x(0)*t1;
      real_type t50 = t8+t11+t47+t49;
      real_type t51 = x(0)*x(3);
      real_type t53 = 0.1214*t51*x(4);
      real_type t59 = x(5)*(0.874E-1*x(0)+0.274436*x(1)+0.874E-1*x(3))*x(2);
      real_type t65 = t47+t49+t13+t15;
      real_type t67 = x(0)+0.157E1*x(1)+x(3);
      real_type t70 = 0.874E-1*x(1)*t67*x(5);
      real_type t76 = x(0)+x(1)+x(2);
      real_type t79 = 0.1214*x(0)*t76*x(4);
      jac[caddr(0,0)] += 0.1214*t2+0.408E-1*x(3);
      jac[caddr(0,1)] += t12;
      jac[caddr(0,2)] += t16;
      jac[caddr(0,3)] += (0.1214*x(0)+0.607E-1*x(1)+0.607E-1*x(2))*t1+0.408E-1*x(0)+0.204E-1*x(1)+0.204E-1*x(2);
      jac[caddr(0,4)] += t31;
      jac[caddr(0,5)] += t34;
      jac[caddr(1,0)] += t12;
      jac[caddr(1,1)] += 0.137218*t10+0.58718E-1*x(2);
      jac[caddr(1,2)] += (0.437E-1*x(0)+0.137218*x(1)+0.437E-1*x(3))*t9+0.187E-1*x(0)+0.58718E-1*x(1)+0.187E-1*x(3);
      jac[caddr(1,3)] += t50;
      jac[caddr(1,4)] += t53;
      jac[caddr(1,5)] += t59;
      jac[caddr(2,0)] += t16;
      jac[caddr(2,1)] += (0.187E-1+0.437E-1*t9)*(x(0)+0.314E1*x(1)+x(3));
      jac[caddr(2,3)] += t65;
      jac[caddr(2,4)] += t53;
      jac[caddr(2,5)] += t70;
      jac[caddr(3,0)] += (0.204E-1+0.607E-1*t1)*(2.0*x(0)+x(1)+x(2));
      jac[caddr(3,1)] += t50;
      jac[caddr(3,2)] += t65;
      jac[caddr(3,4)] += t79;
      jac[caddr(3,5)] += t34;
      jac[caddr(4,0)] += t31;
      jac[caddr(4,1)] += t53;
      jac[caddr(4,2)] += t53;
      jac[caddr(4,3)] += t79;
      jac[caddr(4,4)] += 0.1214*t51*t76;
      jac[caddr(5,0)] += t34;
      jac[caddr(5,1)] += t59;
      jac[caddr(5,2)] += t70;
      jac[caddr(5,3)] += t34;
      jac[caddr(5,5)] += 0.874E-1*t32*t67;
    }
  }

  void
  getExactSolution( dvec_t & x, integer ) const override {
  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getInitialPoint( dvec_t & x, integer ) const override {
    x.fill(1);
    x(0) = 5.54;
    x(1) = 4.4;
    x(2) = 12.02;
    x(3) = 11.82;
    x(4) = 0.702;
    x(5) = 0.852;
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

#define HAS111_BIBTEX \
"@article{Shacham:1990,\n" \
"  author  = {Shacham, Orit and Schacham, Mordechai},\n" \
"  title   = {Finding Boundaries of the Domain of Definition\n" \
"             for Functions Along a One-dimensional Ray},\n" \
"  journal = {ACM Trans. Math. Softw.},\n" \
"  year    = {1990},\n" \
"  volume  = {16},\n" \
"  number  = {3},\n" \
"  pages   = {258--268},\n" \
"  doi     = {10.1145/79505.79511},\n" \
"}\n\n" \
"@book{himmelblau:1972,\n" \
"  author    = {Himmelblau, D.M.},\n" \
"  title     = {Applied nonlinear programming},\n" \
"  year      = {1972},\n" \
"  publisher = {McGraw-Hill}\n" \
"}\n"

class HAS111 : public nonlinearSystem {
  real_type c[10];
  real_type const fact;
public:

  HAS111()
  : nonlinearSystem( "HAS111 function", HAS111_BIBTEX, 10 )
  , fact(1e4)
  {
    c[0] = -6.089;
    c[1] = -17.164;
    c[2] = -34.054;
    c[3] = -5.914;
    c[4] = -24.721;
    c[5] = -14.986;
    c[6] = -24.100;
    c[7] = -10.708;
    c[8] = -26.662;
    c[9] = -22.179;
  }

  real_type
  evalFk( dvec_t const & x, integer k ) const override {
    dvec_t f(n);
    evalF( x, f );
    return f(k);
  }

  void
  evalF( dvec_t const & x, dvec_t & f ) const override {
    real_type ex0 = exp(x(0));
    real_type ex1 = exp(x(1));
    real_type ex2 = exp(x(2));
    real_type ex3 = exp(x(3));
    real_type ex4 = exp(x(4));
    real_type ex5 = exp(x(5));
    real_type ex6 = exp(x(6));
    real_type ex7 = exp(x(7));
    real_type ex8 = exp(x(8));
    real_type ex9 = exp(x(9));
    
    f(0) = 2.0*(ex0+2.0*ex1+2.0*ex2+ex5+ex9-2.0)*ex0;
    f(1) = 4.0*(ex0+2.0*ex1+2.0*ex2+ex5+ex9-2.0)*ex1;
    f(2) = 2.0*ex2*(2.0*ex0+4.0*ex1+5.0*ex2+2.0*ex5+3.0*ex9-5.0+ex6+ex7+2.0*ex8);
    f(3) = 2.0*(ex3+2.0*ex4+ex5+ex6-1.0)*ex3;
    f(4) = 4.0*(ex3+2.0*ex4+ex5+ex6-1.0)*ex4;
    f(5) = 2.0*ex5*(ex0+2.0*ex1+2.0*ex2+2.0*ex5+ex9-3.0+ex3+2.0*ex4+ex6);
    f(6) = 2.0*ex6*(ex3+2.0*ex4+ex5+2.0*ex6-2.0+ex2+ex7+2.0*ex8+ex9);
    f(7) = 2.0*(ex2+ex6+ex7+2.0*ex8+ex9-1.0)*ex7;
    f(8) = 4.0*(ex2+ex6+ex7+2.0*ex8+ex9-1.0)*ex8;
    f(9) = 2.0*ex9*(ex0+2.0*ex1+3.0*ex2+ex5+2.0*ex9-3.0+ex6+ex7+2.0*ex8);

    real_type ss    = ex0+ex1+ex2+ex3+ex4+ex5+ex6+ex7+ex8+ex9;
    real_type logss = log(ss);

    f(0) = ex0*(c[0]+x(0)-logss) + fact*f(0);
    f(1) = ex1*(c[1]+x(1)-logss) + fact*f(1);
    f(2) = ex2*(c[2]+x(2)-logss) + fact*f(2);
    f(3) = ex3*(c[3]+x(3)-logss) + fact*f(3);
    f(4) = ex4*(c[4]+x(4)-logss) + fact*f(4);
    f(5) = ex5*(c[5]+x(5)-logss) + fact*f(5);
    f(6) = ex6*(c[6]+x(6)-logss) + fact*f(6);
    f(7) = ex7*(c[7]+x(7)-logss) + fact*f(7);
    f(8) = ex8*(c[8]+x(8)-logss) + fact*f(8);
    f(9) = ex9*(c[9]+x(9)-logss) + fact*f(9);
  }

  integer
  jacobianNnz() const override
  { return n*n; }

  void
  jacobianPattern( ivec_t & ii, ivec_t & jj ) const override {
    integer kk = 0; // fortran addressing
    for ( integer j = 0; j < n; ++j )
      for ( integer i = 0; i < n; ++i )
        { ii(kk) = i; jj(kk) = j; ++kk; }
  }

  void
  jacobian( dvec_t const & x, dvec_t & jac ) const override {
    real_type ex0 = exp(x(0));
    real_type ex1 = exp(x(1));
    real_type ex2 = exp(x(2));
    real_type ex3 = exp(x(3));
    real_type ex4 = exp(x(4));
    real_type ex5 = exp(x(5));
    real_type ex6 = exp(x(6));
    real_type ex7 = exp(x(7));
    real_type ex8 = exp(x(8));
    real_type ex9 = exp(x(9));

    real_type SS    = ex0+ex1+ex2+ex3+ex4+ex5+ex6+ex7+ex8+ex9;
    real_type logSS = log(SS);

    jac[caddr(0,0)] = ex0*(c[0]+x(0)-logSS)+ex0*(1-ex0/SS);
    jac[caddr(0,1)] = -ex1*ex0/SS;
    jac[caddr(0,2)] = -ex2*ex0/SS;
    jac[caddr(0,3)] = -ex3*ex0/SS;
    jac[caddr(0,4)] = -ex4*ex0/SS;
    jac[caddr(0,5)] = -ex5*ex0/SS;
    jac[caddr(0,6)] = -ex6*ex0/SS;
    jac[caddr(0,7)] = -ex7*ex0/SS;
    jac[caddr(0,8)] = -ex8*ex0/SS;
    jac[caddr(0,9)] = -ex9*ex0/SS;

    jac[caddr(1,0)] = -ex1*ex0/SS;
    jac[caddr(1,1)] = ex1*(c[1]+x(1)-logSS)+ex1*(1-ex1/SS);
    jac[caddr(1,2)] = -ex2*ex1/SS;
    jac[caddr(1,3)] = -ex3*ex1/SS;
    jac[caddr(1,4)] = -ex4*ex1/SS;
    jac[caddr(1,5)] = -ex5*ex1/SS;
    jac[caddr(1,6)] = -ex6*ex1/SS;
    jac[caddr(1,7)] = -ex7*ex1/SS;
    jac[caddr(1,8)] = -ex8*ex1/SS;
    jac[caddr(1,9)] = -ex9*ex1/SS;

    jac[caddr(2,0)] = -ex2*ex0/SS;
    jac[caddr(2,1)] = -ex2*ex1/SS;
    jac[caddr(2,2)] = ex2*(c[2]+x(2)-logSS)+ex2*(1-ex2/SS);
    jac[caddr(2,3)] = -ex3*ex2/SS;
    jac[caddr(2,4)] = -ex4*ex2/SS;
    jac[caddr(2,5)] = -ex5*ex2/SS;
    jac[caddr(2,6)] = -ex6*ex2/SS;
    jac[caddr(2,7)] = -ex7*ex2/SS;
    jac[caddr(2,8)] = -ex8*ex2/SS;
    jac[caddr(2,9)] = -ex9*ex2/SS;

    jac[caddr(3,0)] = -ex3*ex0/SS;
    jac[caddr(3,1)] = -ex3*ex1/SS;
    jac[caddr(3,2)] = -ex3*ex2/SS;
    jac[caddr(3,3)] = ex3*(c[3]+x(3)-logSS)+ex3*(1-ex3/SS);
    jac[caddr(3,4)] = -ex4*ex3/SS;
    jac[caddr(3,5)] = -ex5*ex3/SS;
    jac[caddr(3,6)] = -ex6*ex3/SS;
    jac[caddr(3,7)] = -ex7*ex3/SS;
    jac[caddr(3,8)] = -ex8*ex3/SS;
    jac[caddr(3,9)] = -ex9*ex3/SS;

    jac[caddr(4,0)] = -ex4*ex0/SS;
    jac[caddr(4,1)] = -ex4*ex1/SS;
    jac[caddr(4,2)] = -ex4*ex2/SS;
    jac[caddr(4,3)] = -ex4*ex3/SS;
    jac[caddr(4,4)] = ex4*(c[4]+x(4)-logSS)+ex4*(1-ex4/SS);
    jac[caddr(4,5)] = -ex5*ex4/SS;
    jac[caddr(4,6)] = -ex6*ex4/SS;
    jac[caddr(4,7)] = -ex7*ex4/SS;
    jac[caddr(4,8)] = -ex8*ex4/SS;
    jac[caddr(4,9)] = -ex9*ex4/SS;

    jac[caddr(5,0)] = -ex5*ex0/SS;
    jac[caddr(5,1)] = -ex5*ex1/SS;
    jac[caddr(5,2)] = -ex5*ex2/SS;
    jac[caddr(5,3)] = -ex5*ex3/SS;
    jac[caddr(5,4)] = -ex5*ex4/SS;
    jac[caddr(5,5)] = ex5*(c[5]+x(5)-logSS)+ex5*(1-ex5/SS);
    jac[caddr(5,6)] = -ex6*ex5/SS;
    jac[caddr(5,7)] = -ex7*ex5/SS;
    jac[caddr(5,8)] = -ex8*ex5/SS;
    jac[caddr(5,9)] = -ex9*ex5/SS;

    jac[caddr(6,0)] = -ex6*ex0/SS;
    jac[caddr(6,1)] = -ex6*ex1/SS;
    jac[caddr(6,2)] = -ex6*ex2/SS;
    jac[caddr(6,3)] = -ex6*ex3/SS;
    jac[caddr(6,4)] = -ex6*ex4/SS;
    jac[caddr(6,5)] = -ex6*ex5/SS;
    jac[caddr(6,6)] = ex6*(c[6]+x(6)-logSS)+ex6*(1-ex6/SS);
    jac[caddr(6,7)] = -ex7*ex6/SS;
    jac[caddr(6,8)] = -ex8*ex6/SS;
    jac[caddr(6,9)] = -ex9*ex6/SS;

    jac[caddr(7,0)] = -ex7*ex0/SS;
    jac[caddr(7,1)] = -ex7*ex1/SS;
    jac[caddr(7,2)] = -ex7*ex2/SS;
    jac[caddr(7,3)] = -ex7*ex3/SS;
    jac[caddr(7,4)] = -ex7*ex4/SS;
    jac[caddr(7,5)] = -ex7*ex5/SS;
    jac[caddr(7,6)] = -ex7*ex6/SS;
    jac[caddr(7,7)] = ex7*(c[7]+x(7)-logSS)+ex7*(1-ex7/SS);
    jac[caddr(7,8)] = -ex7*ex8/SS;
    jac[caddr(7,9)] = -ex7*ex9/SS;

    jac[caddr(8,0)] = -ex8*ex0/SS;
    jac[caddr(8,1)] = -ex8*ex1/SS;
    jac[caddr(8,2)] = -ex8*ex2/SS;
    jac[caddr(8,3)] = -ex8*ex3/SS;
    jac[caddr(8,4)] = -ex8*ex4/SS;
    jac[caddr(8,5)] = -ex8*ex5/SS;
    jac[caddr(8,6)] = -ex8*ex6/SS;
    jac[caddr(8,7)] = -ex7*ex8/SS;
    jac[caddr(8,8)] = ex8*(c[8]+x(8)-logSS)+ex8*(1-ex8/SS);
    jac[caddr(8,9)] = -ex8*ex9/SS;

    jac[caddr(9,0)] = -ex9*ex0/SS;
    jac[caddr(9,1)] = -ex9*ex1/SS;
    jac[caddr(9,2)] = -ex9*ex2/SS;
    jac[caddr(9,3)] = -ex9*ex3/SS;
    jac[caddr(9,4)] = -ex9*ex4/SS;
    jac[caddr(9,5)] = -ex9*ex5/SS;
    jac[caddr(9,6)] = -ex9*ex6/SS;
    jac[caddr(9,7)] = -ex7*ex9/SS;
    jac[caddr(9,8)] = -ex8*ex9/SS;
    jac[caddr(9,9)] = ex9*(c[9]+x(9)-logSS)+ex9*(1-ex9/SS);

    jac[caddr(0,0)] += fact*( 2.0*ex0*(2.0*ex0+2.0*ex1+2.0*ex2+ex5+ex9-2.0) );
    jac[caddr(0,1)] += fact*( 4.0*ex1*ex0 );
    jac[caddr(0,2)] += fact*( 4.0*ex2*ex0 );
    jac[caddr(0,5)] += fact*( 2.0*ex5*ex0 );
    jac[caddr(0,9)] += fact*( 2.0*ex9*ex0 );

    jac[caddr(1,0)] += fact*( 4.0*ex1*ex0 );
    jac[caddr(1,1)] += fact*( 4.0*ex1*(4.0*ex1+ex0+2.0*ex2+ex5+ex9-2.0) );
    jac[caddr(1,2)] += fact*( 8.0*ex2*ex1 );
    jac[caddr(1,5)] += fact*( 4.0*ex5*ex1 );
    jac[caddr(1,9)] += fact*( 4.0*ex9*ex1 );

    jac[caddr(2,0)] += fact*( 4.0*ex2*ex0 );
    jac[caddr(2,1)] += fact*( 8.0*ex2*ex1 );
    jac[caddr(2,2)] += fact*( 2.0*ex2*(2.0*ex0+4.0*ex1+10.0*ex2+2.0*ex5+3.0*ex9-5.0+ex6+ex7+2.0*ex8) );
    jac[caddr(2,5)] += fact*( 4.0*ex2*ex5 );
    jac[caddr(2,6)] += fact*( 2.0*ex2*ex6 );
    jac[caddr(2,7)] += fact*( 2.0*ex2*ex7 );
    jac[caddr(2,8)] += fact*( 4.0*ex2*ex8 );
    jac[caddr(2,9)] += fact*( 6.0*ex2*ex9 );

    jac[caddr(3,3)] += fact*( 2.0*ex3*(2.0*ex3+2.0*ex4+ex5+ex6-1.0) );
    jac[caddr(3,4)] += fact*( 4.0*ex4*ex3 );
    jac[caddr(3,5)] += fact*( 2.0*ex5*ex3 );
    jac[caddr(3,6)] += fact*( 2.0*ex6*ex3 );

    jac[caddr(4,3)] += fact*( 4.0*ex4*ex3 );
    jac[caddr(4,4)] += fact*( 4.0*ex4*(4.0*ex4+ex3+ex5+ex6-1.0) );
    jac[caddr(4,5)] += fact*( 4.0*ex5*ex4 );
    jac[caddr(4,6)] += fact*( 4.0*ex6*ex4 );

    jac[caddr(5,0)] += fact*( 2.0*ex5*ex0 );
    jac[caddr(5,1)] += fact*( 4.0*ex5*ex1 );
    jac[caddr(5,2)] += fact*( 4.0*ex2*ex5 );
    jac[caddr(5,3)] += fact*( 2.0*ex5*ex3 );
    jac[caddr(5,4)] += fact*( 4.0*ex5*ex4 );
    jac[caddr(5,5)] += fact*( 2.0*ex5*(ex0+2.0*ex1+2.0*ex2+4.0*ex5+ex9-3.0+ex3+2.0*ex4+ex6) );
    jac[caddr(5,6)] += fact*( 2.0*ex5*ex6 );
    jac[caddr(5,9)] += fact*( 2.0*ex5*ex9 );

    jac[caddr(6,2)] += fact*( 2.0*ex2*ex6 );
    jac[caddr(6,3)] += fact*( 2.0*ex6*ex3 );
    jac[caddr(6,4)] += fact*( 4.0*ex6*ex4 );
    jac[caddr(6,5)] += fact*( 2.0*ex5*ex6 );
    jac[caddr(6,6)] += fact*( 2.0*ex6*(ex3+2.0*ex4+ex5+4.0*ex6-2.0+ex2+ex7+2.0*ex8+ex9) );
    jac[caddr(6,7)] += fact*( 2.0*ex6*ex7 );
    jac[caddr(6,8)] += fact*( 4.0*ex6*ex8 );
    jac[caddr(6,9)] += fact*( 2.0*ex6*ex9 );

    jac[caddr(7,2)] += fact*( 2.0*ex2*ex7 );
    jac[caddr(7,6)] += fact*( 2.0*ex6*ex7 );
    jac[caddr(7,7)] += fact*( 2.0*ex7*(2.0*ex7+ex2+ex6+2.0*ex8+ex9-1.0) );
    jac[caddr(7,8)] += fact*( 4.0*ex7*ex8 );
    jac[caddr(7,9)] += fact*( 2.0*ex7*ex9 );

    jac[caddr(8,2)] += fact*( 4.0*ex2*ex8 );
    jac[caddr(8,6)] += fact*( 4.0*ex6*ex8 );
    jac[caddr(8,7)] += fact*( 4.0*ex7*ex8 );
    jac[caddr(8,8)] += fact*( 4.0*ex8*(4.0*ex8+ex2+ex6+ex7+ex9-1.0) );
    jac[caddr(8,9)] += fact*( 4.0*ex8*ex9 );

    jac[caddr(9,0)] += fact*( 2.0*ex9*ex0 );
    jac[caddr(9,1)] += fact*( 4.0*ex9*ex1 );
    jac[caddr(9,2)] += fact*( 6.0*ex2*ex9 );
    jac[caddr(9,5)] += fact*( 2.0*ex5*ex9 );
    jac[caddr(9,6)] += fact*( 2.0*ex6*ex9 );
    jac[caddr(9,7)] += fact*( 2.0*ex7*ex9 );
    jac[caddr(9,8)] += fact*( 4.0*ex8*ex9 );
    jac[caddr(9,9)] += fact*( 2.0*ex9*(ex0+2.0*ex1+3.0*ex2+ex5+4.0*ex9-3.0+ex6+ex7+2.0*ex8) );

  }

  integer
  numExactSolution() const override
  { return 0; }

  void
  getExactSolution( dvec_t & x, integer ) const override
  { }

  void
  getInitialPoint( dvec_t & x, integer idx ) const override {
    switch ( idx ) {
    case 0:
      x.fill( -1.5-0.5 ); // funziona
      break;
    case 1:
      x.fill( -0.5 ); // no
      break;
    }
  }

  integer
  numInitialPoint() const override
  { return 2; }

  void
  checkIfAdmissible( dvec_t const & x ) const override {
    // funziona solo con questo limite
    //for (  i = 0; i < n; ++i )
    //ASSERT( x(i) < 0 && x(i) > -1, "x range" );
  }

};
