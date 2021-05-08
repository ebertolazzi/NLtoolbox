#include "testsNonlin.hh"
#include <sstream>

namespace NLproblem {

  typedef struct { int_type lo, hi; } stack_node;

  static
  inline
  bool
  GT( int_type i, int_type j, int_type i1, int_type j1 )
  //{ return j==j1 ? i>i1 : j>j1; }
  { return i==i1 ? j>j1 : i>i1; }

  static
  void
  QuickSortIJ(
    ivec_t & I,
    ivec_t & J,
    dvec_t & A,
    int_type total_elems,
    int_type MAX_THRESH = 4
  ) {

    using ::std::swap;

    if ( total_elems <= 1 ) return;
    if ( total_elems <= MAX_THRESH ) goto insert_sort;

    {
      int_type lo = 0;
      int_type hi = total_elems - 1;

      stack_node stack[128];
      stack_node *top = stack;

      while ( top >= stack ) { // Stack not empty

        /* Select median value from among LO, MID, and HI. Rearrange *\
         * LO and HI so the three values are sorted. This lowers the *
         * probability of picking a pathological pivot value and     *
        \* skips a comparison for both the LEFT_PTR and RIGHT_PTR.   */

        int_type & I_hi  = I(hi);
        int_type & I_lo  = I(lo);
        int_type & I_mid = I( (hi + lo) / 2 );

        int_type & J_hi  = J(hi);
        int_type & J_lo  = J(lo);
        int_type & J_mid = J( (hi + lo) / 2 );

        real_type & A_hi  = A(hi);
        real_type & A_lo  = A(lo);
        real_type & A_mid = A( (hi + lo) / 2 );

        if ( GT(I_lo, J_lo, I_mid, J_mid) ) {
          swap(I_mid, I_lo);
          swap(J_mid, J_lo);
          swap(A_mid, A_lo);
        }
        if ( GT(I_mid, J_mid, I_hi, J_hi) ) {
          swap(I_hi, I_mid);
          swap(J_hi, J_mid);
          swap(A_hi, A_mid);
          if ( GT(I_lo, J_lo, I_mid, J_mid) ) {
            swap(I_mid, I_lo);
            swap(J_mid, J_lo);
            swap(A_mid, A_lo);
          }
        }

        int_type IPivot    = I_mid;
        int_type JPivot    = J_mid;
        int_type left_ptr  = lo + 1;
        int_type right_ptr = hi - 1;

        /* Here's the famous ``collapse the walls'' section of quicksort. *\
         * Gotta like those tight inner loops!  They are the main reason  *
        \* that this algorithm runs much faster than others.              */

        do {

          while ( GT(IPivot, JPivot, I(left_ptr), J(left_ptr) )  ) ++left_ptr;
          while ( GT(I(right_ptr), J(right_ptr), IPivot, JPivot) ) --right_ptr;

          if ( left_ptr < right_ptr ) {
            swap( I(left_ptr), I(right_ptr) );
            swap( J(left_ptr), J(right_ptr) );
            swap( A(left_ptr), A(right_ptr) );
            ++left_ptr;
            --right_ptr;
          } else if ( left_ptr == right_ptr ) {
            ++left_ptr;
            --right_ptr;
            break;
          }

        } while ( left_ptr <= right_ptr );

        /* Set up pointers for next iteration.  First determine whether   *\
         * left and right partitions are below the threshold size. If so, *
         * ignore one or both.  Otherwise, push the larger partition's    *
        \* bounds on the stack and continue sorting the smaller one.      */

        if ( (right_ptr - lo) <= MAX_THRESH ) {
          if ((hi - left_ptr) <= MAX_THRESH ) {
            --top;
            lo = top -> lo;
            hi = top -> hi;
          } else {
            lo = left_ptr;
          }
        } else if ((hi - left_ptr) <= MAX_THRESH) {
          hi = right_ptr;
        } else if ((right_ptr - lo) > (hi - left_ptr)) {
          top -> lo = lo;
          top -> hi = right_ptr;
          ++top;
          lo = left_ptr;
        } else {
          top -> lo = left_ptr;
          top -> hi = hi;
          ++top;
          hi = right_ptr;
        }
      }

      /* Once the BASE_PTR array is partially sorted by quicksort the rest  *\
       * is completely sorted using insertion sort, since this is efficient *
       * for partitions below MAX_THRESH size.                              */
    }

  insert_sort:

    for ( int_type i = 1; i < total_elems; ++i ) {
      for ( int_type j = i; j > 0; --j ) {
        if ( GT( I(j), J(j), I(j-1), J(j-1) ) ) break;
        swap( I(j), I(j-1) );
        swap( J(j), J(j-1) );
        swap( A(j), A(j-1) );
      }
    }
  }

  int_type
  nonlinearSystem::fill_CSR(
    dvec_t const & x,
    ivec_t       & R,
    ivec_t       & J,
    dvec_t       & values
  ) const {
    int_type const nnz = jacobianNnz();
    ivec_t I( nnz );
    J.resize( nnz );
    values.resize( nnz );
    jacobianPattern( I, J );
    jacobian( x, values );
    // reorder
    QuickSortIJ( I, J, values, nnz );
    R.resize( numEqns()+1 );
    R.setZero();
    for ( int_type i = 0; i < nnz;       ++i ) ++R(I(i)+1);
    for ( int_type i = 1; i < numEqns(); ++i ) R(i+1) += R(i);
    return nnz;
  };

  #include "tests/ArtificialTestOfNowakAndWeimann.cxx"
  #include "tests/BadlyScaledAugmentedPowellFunction.cxx"
  #include "tests/Beale.cxx"
  #include "tests/Bertolazzi.cxx"
  #include "tests/BiggsEXPfunctions.cxx"
  #include "tests/BoggsFunction.cxx"
  #include "tests/Bohachevsky.cxx"
  #include "tests/BoxAndBettsExponentialQuadraticSum.cxx"
  #include "tests/BoxProblem.cxx"
  #include "tests/Box3.cxx"
  #include "tests/BraninRCOS.cxx"
  #include "tests/BrownAlmostLinearFunction.cxx"
  #include "tests/BrownAndConteFunction.cxx"
  #include "tests/BrownAndDennis.cxx"
  #include "tests/BrownAndGearhartFunction.cxx"
  #include "tests/BrownBadlyScaled.cxx"
  #include "tests/BrownFunction.cxx"
  #include "tests/BroydenBandedFunction.cxx"
  #include "tests/BroydenTridiagonalFunction.cxx"
  #include "tests/BUNLSI.cxx"
  #include "tests/BurdenAndFaires.cxx"
  #include "tests/Chandrasekhar.cxx"
  #include "tests/ChebyquadFunction.cxx"
  #include "tests/ChemicalEquilibriumApplication.cxx"
  #include "tests/CliffFunction.cxx"
  #include "tests/Colville.cxx"
  #include "tests/CombustionApplication.cxx"
  #include "tests/ComplementaryFunction.cxx"
  #include "tests/CompressibilityFactorFromTheRKequation.cxx"
  #include "tests/CountercurrentReactorsProblem.cxx"
  #include "tests/CraggAndLevyProblem.cxx"
  #include "tests/CubeFunction.cxx"
  #include "tests/DarvishiBarati.cxx"
  #include "tests/DeistSefor.cxx"
  #include "tests/DennisAndGay.cxx"
  #include "tests/DennisAndSchnabel2x2example.cxx"
  #include "tests/DeVilliersGlasser.cxx"
  #include "tests/DiagonalFunctionMulQO.cxx"
  #include "tests/DiscreteBoundaryValueFunction.cxx"
  #include "tests/DiscreteIntegralEquationFunction.cxx"
  #include "tests/DixonFunction.cxx"
  #include "tests/Easom.cxx"
  #include "tests/EsterificReaction.cxx"
  #include "tests/ExponentialFunction.cxx"
  #include "tests/ExponentialSine.cxx"
  #include "tests/ExtendedEigerSikorskiStenger.cxx"
  #include "tests/ExtendedKearfottFunction.cxx"
  #include "tests/ExtendedPowellSingularFunction.cxx"
  #include "tests/FreudensteinRothFunction.cxx"
  #include "tests/Function15.cxx"
  #include "tests/Function18.cxx"
  #include "tests/Function21.cxx"
  #include "tests/Function27.cxx"
  #include "tests/Gauss.cxx"
  #include "tests/GeneralizedRosenbrock.cxx"
  #include "tests/GeometricProgrammingFunction.cxx"
  #include "tests/GheriMancino.cxx"
  #include "tests/GoldsteinPrice.cxx"
  #include "tests/GregoryAndKarney.cxx"
  #include "tests/GriewankFunction.cxx"
  #include "tests/Gulf.cxx"
  #include "tests/Hammarling.cxx"
  #include "tests/HanbookFunction.cxx"
  #include "tests/HanSunHan.cxx"
  #include "tests/HAS.cxx"
  #include "tests/HelicalValleyFunction.cxx"
  #include "tests/HiebertChem.cxx"
  #include "tests/Hilbert.cxx"
  #include "tests/Himmelblau.cxx"
  #include "tests/InfRefluxFunction.cxx"
  #include "tests/IntervalArithmeticBenchmarks.cxx"
  #include "tests/JennrichAndSampsonFunction.cxx"
  #include "tests/KelleyFunction.cxx"
  #include "tests/KinematicApplication.cxx"
  #include "tests/Leon.cxx"
  #include "tests/LinearFunctionFullRank.cxx"
  #include "tests/LinearFunctionRank1.cxx"
  #include "tests/LogarithmicFunction.cxx"
  #include "tests/McCormicFunction.cxx"
  #include "tests/McKinnon.cxx"
  #include "tests/MexicanHatFunction.cxx"
  #include "tests/MieleAndCantrellFunction.cxx"
  #include "tests/NonlinearIntegralEquations.cxx"
  #include "tests/Order10to11function.cxx"
  #include "tests/PavianiFunction.cxx"
  #include "tests/Penalty.cxx"
  #include "tests/PowellBadlyScaledFunction.cxx"
  #include "tests/PowellQuarticFunction.cxx"
  #include "tests/Powell3D.cxx"
  #include "tests/RooseKullaLombMeressoo.cxx"
  #include "tests/SampleProblem18.cxx"
  #include "tests/SampleProblem19.cxx"
  #include "tests/ScalarProblem.cxx"
  #include "tests/Schaffer.cxx"
  #include "tests/SchubertBroydenFunction.cxx"
  #include "tests/Semiconductor2D.cxx"
  #include "tests/Shacham.cxx"
  #include "tests/Shekel.cxx"
  #include "tests/ShenYpma.cxx"
  #include "tests/Shubert.cxx"
  #include "tests/SingularFunction.cxx"
  #include "tests/SingularSystem.cxx"
  #include "tests/SIRtest.cxx"
  #include "tests/SixHumpCamelBackFunction.cxx"
  #include "tests/SoniaKrzyworzcka.cxx"
  #include "tests/Spedicato.cxx"
  #include "tests/SSTnonlinearityTerm.cxx"
  #include "tests/StrictlyConvexFunction.cxx"
  #include "tests/Toint.cxx"
  #include "tests/TridimensionalValley.cxx"
  #include "tests/TrigonometricFunction.cxx"
  #include "tests/TrigonometricExponentialSystem.cxx"
  #include "tests/TroeschFunction.cxx"
  #include "tests/TwoPointBoundaryValueProblem.cxx"
  #include "tests/VariablyDimensionedFunction.cxx"
  #include "tests/Weibull.cxx"
  #include "tests/WoodFunction.cxx"
  #include "tests/WatsonFunction.cxx"
  #include "tests/XiaoYin.cxx"
  #include "tests/YixunShi.cxx"
  #include "tests/ZeroJacobianFunction.cxx"

  std::vector<nonlinearSystem*> theProblems;
  std::map<string,int_type>     theProblemsMap;

  void
  initProblems() {

    theProblems.push_back( new ArtificialTestOfNowakAndWeimann() );
    theProblems.push_back( new BadlyScaledAugmentedPowellFunction(3) );
    theProblems.push_back( new BadlyScaledAugmentedPowellFunction(30) );
    theProblems.push_back( new BadlyScaledAugmentedPowellFunction(300) );
    theProblems.push_back( new BadlyScaledAugmentedPowellFunction(3000) );
    theProblems.push_back( new Beale() );
    theProblems.push_back( new BertolazziRootPlusSquare() );
    theProblems.push_back( new BertolazziAtanPlusQuadratic() );
    theProblems.push_back( new BertolazziHard() );
    theProblems.push_back( new BertolazziSingleEQ() );

    theProblems.push_back( new BiggsEXP2function() );
    theProblems.push_back( new BiggsEXP3function() );
    theProblems.push_back( new BiggsEXP4function() );
    theProblems.push_back( new BiggsEXP5function() );
    theProblems.push_back( new BiggsEXP6function() );
    theProblems.push_back( new BoggsFunction() );
    theProblems.push_back( new BohachevskyN1() );
    theProblems.push_back( new BohachevskyN2() );
    theProblems.push_back( new BohachevskyN3() );

    theProblems.push_back( new BoxAndBettsExponentialQuadraticSum() );
    theProblems.push_back( new BoxProblem() );
    theProblems.push_back( new Box3() );
    theProblems.push_back( new BraninRCOS() );
    theProblems.push_back( new BrownAlmostLinearFunction(5) );
    theProblems.push_back( new BrownAlmostLinearFunction(15) );
    theProblems.push_back( new BrownAlmostLinearFunction(25) );
    theProblems.push_back( new BrownAndConteFunction() );
    theProblems.push_back( new BrownAndDennis() );
    theProblems.push_back( new BrownAndGearhartFunction() );
    theProblems.push_back( new BrownBadlyScaled() );
    theProblems.push_back( new BrownFunction() );
    theProblems.push_back( new BroydenBandedFunction() );

    theProblems.push_back( new BroydenTridiagonalFunction(0.1,1,5) );
    theProblems.push_back( new BroydenTridiagonalFunction(0.1,1,10) );
    theProblems.push_back( new BroydenTridiagonalFunction(0.1,1,500) );
    theProblems.push_back( new BroydenTridiagonalFunction(0.5,1,5) );
    theProblems.push_back( new BroydenTridiagonalFunction(0.5,1,10) );
    theProblems.push_back( new BroydenTridiagonalFunction(0.5,1,500) );

    theProblems.push_back( new BUNLSI5() );
    theProblems.push_back( new BUNLSI6() );

    theProblems.push_back( new BurdenAndFaires() );

    theProblems.push_back( new Chandrasekhar(0.9999,10) );
    theProblems.push_back( new Chandrasekhar(0.9999,50) );
    theProblems.push_back( new Chandrasekhar(0.9,10) );
    theProblems.push_back( new Chandrasekhar(0.9,50) );

    theProblems.push_back( new ChebyquadFunction(1) );
    theProblems.push_back( new ChebyquadFunction(2) );
    theProblems.push_back( new ChebyquadFunction(3) );
    theProblems.push_back( new ChebyquadFunction(4) );
    theProblems.push_back( new ChebyquadFunction(5) );
    theProblems.push_back( new ChebyquadFunction(6) );
    theProblems.push_back( new ChebyquadFunction(7) );
    theProblems.push_back( new ChebyquadFunction(8) );
    theProblems.push_back( new ChebyquadFunction(9) );

    theProblems.push_back( new ChemicalEquilibriumApplication() );
    theProblems.push_back( new ChemicalEquilibriumPartialMethaneOxidation() );
    theProblems.push_back( new ChemicalReactorEquilibriumConversion() );
    theProblems.push_back( new ChemicalReactorSteadyState() );

    theProblems.push_back( new CliffFunction() );
    theProblems.push_back( new Colville() );
    theProblems.push_back( new CombustionApplication() );

    theProblems.push_back( new ComplementaryFunction(2) );
    theProblems.push_back( new ComplementaryFunction(16) );
    theProblems.push_back( new ComplementaryFunction(128) );

    theProblems.push_back( new CompressibilityFactorFromTheRKequation() );

    theProblems.push_back( new CountercurrentReactorsProblem1(6) );
    theProblems.push_back( new CountercurrentReactorsProblem1(12) );
    theProblems.push_back( new CountercurrentReactorsProblem1(50) );
    theProblems.push_back( new CountercurrentReactorsProblem2(6) );
    theProblems.push_back( new CountercurrentReactorsProblem2(12) );
    theProblems.push_back( new CountercurrentReactorsProblem2(50) );
    theProblems.push_back( new CraggAndLevyProblem() );
    theProblems.push_back( new CubeFunction() );

    theProblems.push_back( new CutlipsSteadyStateForReactionRateEquations(0) );
    theProblems.push_back( new CutlipsSteadyStateForReactionRateEquations(1) );
    theProblems.push_back( new CutlipsSteadyStateForReactionRateEquations(2) );

    theProblems.push_back( new DarvishiBarati() );

    theProblems.push_back( new DennisAndGay6eqN1() );
    theProblems.push_back( new DennisAndGay6eqN2() );
    theProblems.push_back( new DennisAndGay6eqN3() );
    theProblems.push_back( new DennisAndGay6eqN4() );
    theProblems.push_back( new DennisAndGay6eqN5() );
    theProblems.push_back( new DennisAndGay8eqN1() );
    theProblems.push_back( new DennisAndGay8eqN2() );
    theProblems.push_back( new DennisAndGay8eqN3() );
    theProblems.push_back( new DennisAndGay8eqN4() );
    theProblems.push_back( new DennisAndGay8eqN5() );

    theProblems.push_back( new DennisAndSchnabel2x2example() );

    theProblems.push_back( new DeVilliersGlasser01() );
    theProblems.push_back( new DeVilliersGlasser02() );

    theProblems.push_back( new DiagonalFunctionMulQO(9) );
    theProblems.push_back( new DiagonalFunctionMulQO(27) );

    theProblems.push_back( new DiscreteBoundaryValueFunction(10) );
    theProblems.push_back( new DiscreteBoundaryValueFunction(50) );
    theProblems.push_back( new DiscreteBoundaryValueFunction(100) );
    theProblems.push_back( new DiscreteBoundaryValueFunction(500) );
    theProblems.push_back( new DiscreteBoundaryValueFunction(5000) );

    theProblems.push_back( new DiscreteIntegralEquationFunction(2) );
    theProblems.push_back( new DiscreteIntegralEquationFunction(5) );
    theProblems.push_back( new DiscreteIntegralEquationFunction(10) );
    theProblems.push_back( new DiscreteIntegralEquationFunction(100) );

    theProblems.push_back( new DixonFunction(80) );
    theProblems.push_back( new DixonFunction(2000) );
    theProblems.push_back( new DixonFunction(5000) );

    theProblems.push_back( new Easom() );
    theProblems.push_back( new EsterificReaction() );

    theProblems.push_back( new ExponentialFunction1(2) );
    theProblems.push_back( new ExponentialFunction1(10) );
    theProblems.push_back( new ExponentialFunction1(50) );
    theProblems.push_back( new ExponentialFunction1(500) );
    theProblems.push_back( new ExponentialFunction1(5000) );
    theProblems.push_back( new ExponentialFunction2(2) );
    theProblems.push_back( new ExponentialFunction2(10) );
    theProblems.push_back( new ExponentialFunction2(50) );
    theProblems.push_back( new ExponentialFunction2(500) );
    theProblems.push_back( new ExponentialFunction2(5000) );
    theProblems.push_back( new ExponentialFunction3(2) );
    theProblems.push_back( new ExponentialFunction3(10) );
    theProblems.push_back( new ExponentialFunction3(50) );
    theProblems.push_back( new ExponentialFunction3(500) );
    theProblems.push_back( new ExponentialFunction3(5000) );

    theProblems.push_back( new ExponentialSine() );
    theProblems.push_back( new ExtendedEigerSikorskiStenger() );
    theProblems.push_back( new ExtendedKearfottFunction() );
    theProblems.push_back( new ExtendedPowellSingularFunction() );

    theProblems.push_back( new FractionalConversionInAchemicalReactor() );
    theProblems.push_back( new FractionalConversionInAchemicalReactor2() );

    theProblems.push_back( new FreudensteinRothFunction() );

    theProblems.push_back( new Function15(10) );
    theProblems.push_back( new Function15(50) );
    theProblems.push_back( new Function15(100) );

    theProblems.push_back( new Function18(3) );
    theProblems.push_back( new Function18(9) );
    theProblems.push_back( new Function18(27) );

    theProblems.push_back( new Function21(3) );
    theProblems.push_back( new Function21(6) );
    theProblems.push_back( new Function21(36) );

    theProblems.push_back( new Function27(2) );
    theProblems.push_back( new Function27(10) );
    theProblems.push_back( new Function27(100) );

    theProblems.push_back( new Gauss() );

    theProblems.push_back( new GeneralizedRosenbrock(2) );
    theProblems.push_back( new GeneralizedRosenbrock(10) );
    theProblems.push_back( new GeneralizedRosenbrock(50) );
    theProblems.push_back( new GeneralizedRosenbrock(500) );

    theProblems.push_back( new GeometricProgrammingFunction(5) );
    theProblems.push_back( new GeometricProgrammingFunction(10) );
    theProblems.push_back( new GeometricProgrammingFunction(50) );
    theProblems.push_back( new GeometricProgrammingFunction(100) );

    theProblems.push_back( new GheriMancino(10) );
    theProblems.push_back( new GheriMancino(30) );
    theProblems.push_back( new GheriMancino(100) );

    theProblems.push_back( new GoldsteinPrice() );
    theProblems.push_back( new GregoryAndKarney(10) );
    theProblems.push_back( new GregoryAndKarney(20) );
    theProblems.push_back( new GregoryAndKarney(100) );

    theProblems.push_back( new GriewankFunction(2) );
    theProblems.push_back( new GriewankFunction(5) );
    theProblems.push_back( new GriewankFunction(10) );

    theProblems.push_back( new Gulf() );

    theProblems.push_back( new Hammarling2x2matrixSquareRoot() );
    theProblems.push_back( new Hammarling3x3matrixSquareRootProblemN1() );
    theProblems.push_back( new Hammarling3x3matrixSquareRootProblemN2() );
    theProblems.push_back( new Hammarling3x3matrixSquareRootProblemN3() );

    theProblems.push_back( new HanbookFunction(2) );
    theProblems.push_back( new HanbookFunction(10) );
    theProblems.push_back( new HanbookFunction(50) );

    theProblems.push_back( new HanSunHan() );

    theProblems.push_back( new HAS64(1e2)  );
    theProblems.push_back( new HAS64(1e4)  );
    theProblems.push_back( new HAS64(1e6)  );
    theProblems.push_back( new HAS64(1e8)  );
    theProblems.push_back( new HAS64(1e10) );
    theProblems.push_back( new HAS93(1e2)  );
    theProblems.push_back( new HAS93(1e4)  );
    theProblems.push_back( new HAS93(1e6)  );
    theProblems.push_back( new HAS93(1e8)  );
    theProblems.push_back( new HAS93(1e10) );
    theProblems.push_back( new HAS111() );

    theProblems.push_back( new HelicalValleyFunction() );

    theProblems.push_back( new HiebertChem10x10() );
    theProblems.push_back( new HiebertChem2x2() );
    theProblems.push_back( new HiebertChem6x6() );

    theProblems.push_back( new Hiebert3ChemicalEquilibriumProblem(10) );
    theProblems.push_back( new Hiebert3ChemicalEquilibriumProblem(40) );

    theProblems.push_back( new Hilbert(4) );
    theProblems.push_back( new Hilbert(8) );
    theProblems.push_back( new Hilbert(32) );
    theProblems.push_back( new Hilbert(64) );
    theProblems.push_back( new Hilbert(256) );

    theProblems.push_back( new Himmelblau() );
    theProblems.push_back( new InfRefluxFunction() );
    theProblems.push_back( new IntervalArithmeticBenchmarks() );
    theProblems.push_back( new JennrichAndSampsonFunction() );
    theProblems.push_back( new KelleyFunction() );
    theProblems.push_back( new KinematicApplication() );
    theProblems.push_back( new Leon() );
    theProblems.push_back( new LinearFunctionFullRank() );
    theProblems.push_back( new LinearFunctionRank1() );

    theProblems.push_back( new LogarithmicFunction(2) );
    theProblems.push_back( new LogarithmicFunction(10) );
    theProblems.push_back( new LogarithmicFunction(50) );
    theProblems.push_back( new LogarithmicFunction(500) );
    theProblems.push_back( new LogarithmicFunction(5000) );

    theProblems.push_back( new McCormicFunction() );
    theProblems.push_back( new McKinnon() );

    theProblems.push_back( new MexicanHatFunction(1e2) );
    theProblems.push_back( new MexicanHatFunction(1e4) );
    theProblems.push_back( new MexicanHatFunction(1e6) );
    theProblems.push_back( new MexicanHatFunction(1e8) );

    theProblems.push_back( new MieleAndCantrellFunction() );

    theProblems.push_back( new ModelEquationsForTheCSTR() );
    theProblems.push_back( new ModelEquationsForCombustionOfPropane1() );
    theProblems.push_back( new ModelEquationsForCombustionOfPropane2() );

    theProblems.push_back( new NonlinearIntegralEquations() );
    theProblems.push_back( new Order10to11function() );

    theProblems.push_back( new PavianiFunction() );

    theProblems.push_back( new PenaltyIfunction(10) );
    theProblems.push_back( new PenaltyIfunction(50) );

    theProblems.push_back( new PenaltyN1(2) );
    theProblems.push_back( new PenaltyN1(10) );
    theProblems.push_back( new PenaltyN1(50) );
    theProblems.push_back( new PenaltyN1(200) );

    theProblems.push_back( new PenaltyN2(2) );
    theProblems.push_back( new PenaltyN2(10) );
    theProblems.push_back( new PenaltyN2(50) );
    theProblems.push_back( new PenaltyN2(200) );

    theProblems.push_back( new PipelineNetworkProblem() );
    theProblems.push_back( new PipelineNetworkProblem2() );

    theProblems.push_back( new PowellBadlyScaledFunction() );
    theProblems.push_back( new PowellQuarticFunction() );
    theProblems.push_back( new Powell3D() );

    theProblems.push_back( new RooseKullaLombMeressoo129() );

    theProblems.push_back( new RooseKullaLombMeressoo201(10) );
    theProblems.push_back( new RooseKullaLombMeressoo201(100) );
    theProblems.push_back( new RooseKullaLombMeressoo201(500) );

    theProblems.push_back( new RooseKullaLombMeressoo202(10) );
    theProblems.push_back( new RooseKullaLombMeressoo202(100) );

    theProblems.push_back( new RooseKullaLombMeressoo203(2) );
    theProblems.push_back( new RooseKullaLombMeressoo203(5) );
    theProblems.push_back( new RooseKullaLombMeressoo203(7) );
    theProblems.push_back( new RooseKullaLombMeressoo203(10) );
    theProblems.push_back( new RooseKullaLombMeressoo203(15) );
    theProblems.push_back( new RooseKullaLombMeressoo203(20) );
    theProblems.push_back( new RooseKullaLombMeressoo203(30) );
    theProblems.push_back( new RooseKullaLombMeressoo203(50) );

    theProblems.push_back( new RooseKullaLombMeressoo204(10) );
    theProblems.push_back( new RooseKullaLombMeressoo204(100) );
    theProblems.push_back( new RooseKullaLombMeressoo204(500) );

    theProblems.push_back( new RooseKullaLombMeressoo205(2) );
    theProblems.push_back( new RooseKullaLombMeressoo205(5) );
    theProblems.push_back( new RooseKullaLombMeressoo205(10) );
    theProblems.push_back( new RooseKullaLombMeressoo205(20) );
    theProblems.push_back( new RooseKullaLombMeressoo205(30) );

    theProblems.push_back( new RooseKullaLombMeressoo206(2) );
    theProblems.push_back( new RooseKullaLombMeressoo206(5) );
    theProblems.push_back( new RooseKullaLombMeressoo206(10) );
    theProblems.push_back( new RooseKullaLombMeressoo206(20) );
    theProblems.push_back( new RooseKullaLombMeressoo206(30) );

    theProblems.push_back( new RooseKullaLombMeressoo207(2) );
    theProblems.push_back( new RooseKullaLombMeressoo207(5) );
    theProblems.push_back( new RooseKullaLombMeressoo207(10) );
    theProblems.push_back( new RooseKullaLombMeressoo207(20) );
    theProblems.push_back( new RooseKullaLombMeressoo207(30) );

    theProblems.push_back( new RooseKullaLombMeressoo208(2) );
    theProblems.push_back( new RooseKullaLombMeressoo208(5) );
    theProblems.push_back( new RooseKullaLombMeressoo208(10) );
    theProblems.push_back( new RooseKullaLombMeressoo208(20) );
    theProblems.push_back( new RooseKullaLombMeressoo208(30) );

    theProblems.push_back( new RooseKullaLombMeressoo209(10) );
    theProblems.push_back( new RooseKullaLombMeressoo209(100) );
    theProblems.push_back( new RooseKullaLombMeressoo209(500) );

    theProblems.push_back( new RooseKullaLombMeressoo210(10) ); // no solution
    theProblems.push_back( new RooseKullaLombMeressoo211(10) ); // no solution

    theProblems.push_back( new RooseKullaLombMeressoo212(10) );
    theProblems.push_back( new RooseKullaLombMeressoo212(100) );
    theProblems.push_back( new RooseKullaLombMeressoo212(500) );

    theProblems.push_back( new RooseKullaLombMeressoo213(2) );
    theProblems.push_back( new RooseKullaLombMeressoo213(5) );
    theProblems.push_back( new RooseKullaLombMeressoo213(10) );
    theProblems.push_back( new RooseKullaLombMeressoo213(20) );
    theProblems.push_back( new RooseKullaLombMeressoo213(30) );

    theProblems.push_back( new RooseKullaLombMeressoo214(2) );
    theProblems.push_back( new RooseKullaLombMeressoo214(5) );
    theProblems.push_back( new RooseKullaLombMeressoo214(10) );
    theProblems.push_back( new RooseKullaLombMeressoo214(20) );
    theProblems.push_back( new RooseKullaLombMeressoo214(30) );

    theProblems.push_back( new RooseKullaLombMeressoo215(2) );
    theProblems.push_back( new RooseKullaLombMeressoo215(5) );
    theProblems.push_back( new RooseKullaLombMeressoo215(6) );
    theProblems.push_back( new RooseKullaLombMeressoo215(10) );
    theProblems.push_back( new RooseKullaLombMeressoo215(20) );
    theProblems.push_back( new RooseKullaLombMeressoo215(30) );

    theProblems.push_back( new RooseKullaLombMeressoo216(2) );
    theProblems.push_back( new RooseKullaLombMeressoo216(5) );
    theProblems.push_back( new RooseKullaLombMeressoo216(7) );
    theProblems.push_back( new RooseKullaLombMeressoo216(9) );

    theProblems.push_back( new RooseKullaLombMeressoo217(2) );
    theProblems.push_back( new RooseKullaLombMeressoo217(5) );
    theProblems.push_back( new RooseKullaLombMeressoo217(10) );
    theProblems.push_back( new RooseKullaLombMeressoo217(20) );
    theProblems.push_back( new RooseKullaLombMeressoo217(30) );

    theProblems.push_back( new RooseKullaLombMeressoo218(2) );
    theProblems.push_back( new RooseKullaLombMeressoo218(5) );
    theProblems.push_back( new RooseKullaLombMeressoo218(10) );
    theProblems.push_back( new RooseKullaLombMeressoo218(20) );
    theProblems.push_back( new RooseKullaLombMeressoo218(30) );

    theProblems.push_back( new RooseKullaLombMeressoo219(2) );
    theProblems.push_back( new RooseKullaLombMeressoo219(5) );
    theProblems.push_back( new RooseKullaLombMeressoo219(10) );
    theProblems.push_back( new RooseKullaLombMeressoo219(20) );
    theProblems.push_back( new RooseKullaLombMeressoo219(30) );

    theProblems.push_back( new SampleProblem18() );
    theProblems.push_back( new SampleProblem19() );
    theProblems.push_back( new ScalarProblem() );

    theProblems.push_back( new SchafferF6() );
    theProblems.push_back( new SchafferF7() );

    theProblems.push_back( new SchubertBroydenFunction(10) );
    theProblems.push_back( new SchubertBroydenFunction(100) );
    theProblems.push_back( new SchubertBroydenFunction(1000) );
    theProblems.push_back( new SchubertBroydenFunction(5000) );

    theProblems.push_back( new Semiconductor2D() );

    theProblems.push_back( new ShekelSQRN5() );
    theProblems.push_back( new ShekelSQRN7() );
    theProblems.push_back( new ShekelSQRN10() );

    theProblems.push_back( new ShenYpma5() );
    theProblems.push_back( new ShenYpma7() );
    theProblems.push_back( new ShenYpma8() );

    theProblems.push_back( new Shubert() );

    theProblems.push_back( new SingularFunction(2) );
    theProblems.push_back( new SingularFunction(10) );
    theProblems.push_back( new SingularFunction(50) );
    theProblems.push_back( new SingularFunction(500) );
    theProblems.push_back( new SingularFunction(5000) );

    theProblems.push_back( new SingularSystemA() );
    theProblems.push_back( new SingularSystemB() );
    theProblems.push_back( new SingularSystemC() );
    theProblems.push_back( new SingularSystemD() );
    theProblems.push_back( new SingularSystemE() );
    theProblems.push_back( new SingularSystemF() );

    theProblems.push_back( new SingularSystemP2() );
    theProblems.push_back( new SingularSystemP3() );
    theProblems.push_back( new SingularSystemP4() );
    theProblems.push_back( new SingularSystemP5() );
    theProblems.push_back( new SingularSystemP6() );
    theProblems.push_back( new SingularSystemP7() );
    theProblems.push_back( new SingularSystemP8() );
    theProblems.push_back( new SingularSystemP9() );

    theProblems.push_back( new SIRtest(2) );
    theProblems.push_back( new SIRtest(20) );
    theProblems.push_back( new SIRtest(100) );

    theProblems.push_back( new SixHumpCamelBackFunction() );

    theProblems.push_back( new SoniaKrzyworzcka1() );
    theProblems.push_back( new SoniaKrzyworzcka2() );

    theProblems.push_back( new SpedicatoFunction17(10) );
    theProblems.push_back( new SpedicatoFunction17(50) );
    theProblems.push_back( new SpedicatoFunction17(100) );
    theProblems.push_back( new SpedicatoFunction17(500) );

    theProblems.push_back( new SSTnonlinearityTerm(0) );
    theProblems.push_back( new SSTnonlinearityTerm(1) );

    theProblems.push_back( new StrictlyConvexFunction1(2) );
    theProblems.push_back( new StrictlyConvexFunction1(10) );
    theProblems.push_back( new StrictlyConvexFunction1(50) );
    theProblems.push_back( new StrictlyConvexFunction1(500) );
    theProblems.push_back( new StrictlyConvexFunction1(5000) );

    theProblems.push_back( new StrictlyConvexFunction2(2) );
    theProblems.push_back( new StrictlyConvexFunction2(10) );
    theProblems.push_back( new StrictlyConvexFunction2(50) );
    theProblems.push_back( new StrictlyConvexFunction2(500) );
    theProblems.push_back( new StrictlyConvexFunction2(5000) );

    theProblems.push_back( new Toint225(10) );
    theProblems.push_back( new Toint225(100) );
    theProblems.push_back( new Toint225(500) );

    theProblems.push_back( new TridimensionalValley() );

    theProblems.push_back( new TrigonometricExponentialSystem1(10) );
    theProblems.push_back( new TrigonometricExponentialSystem1(50) );
    theProblems.push_back( new TrigonometricExponentialSystem1(500) );

    theProblems.push_back( new TrigonometricExponentialSystem2(9) );
    theProblems.push_back( new TrigonometricExponentialSystem2(27) );
    theProblems.push_back( new TrigonometricExponentialSystem2(81) );

    theProblems.push_back( new TrigExp(10) );
    theProblems.push_back( new TrigExp(100) );

    theProblems.push_back( new TrigonometricFunction(2) );
    theProblems.push_back( new TrigonometricFunction(10) );
    theProblems.push_back( new TrigonometricFunction(50) );

    theProblems.push_back( new TroeschFunction(2) );
    theProblems.push_back( new TroeschFunction(10) );
    theProblems.push_back( new TroeschFunction(50) );
    theProblems.push_back( new TroeschFunction(500) );
    theProblems.push_back( new TroeschFunction(5000) );

    theProblems.push_back( new TwoPointBoundaryValueProblem(10) );
    theProblems.push_back( new TwoPointBoundaryValueProblem(100) );
    theProblems.push_back( new TwoPointBoundaryValueProblem(1000) );
    theProblems.push_back( new TwoPointBoundaryValueProblem(5000) );

    theProblems.push_back( new VariablyDimensionedFunction(5) );
    theProblems.push_back( new VariablyDimensionedFunction(10) );
    theProblems.push_back( new VariablyDimensionedFunction(50) );

    theProblems.push_back( new WatsonFunction() );
    theProblems.push_back( new Weibull() );

    theProblems.push_back( new WoodFunction() );

    theProblems.push_back( new XiaoYin1() );
    theProblems.push_back( new XiaoYin2() );
    theProblems.push_back( new XiaoYin3() );

    theProblems.push_back( new YixunShi1() );
    theProblems.push_back( new YixunShi2() );
    theProblems.push_back( new YixunShi3() );
    theProblems.push_back( new YixunShi4() );

    theProblems.push_back( new ZeroJacobianFunction(10) );
    theProblems.push_back( new ZeroJacobianFunction(50) );
    theProblems.push_back( new ZeroJacobianFunction(101) );

    for ( int_type i = 0; i < int_type( theProblemsMap.size() ); ++i ) {
      theProblemsMap[theProblems[i]->title()] = i;
    }

  }

}
