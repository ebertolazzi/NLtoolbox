#include "testsNonlin.hh"
#include "mex_utils.hh"

#include <map>
#include <vector>

#define MEX_ERROR_MESSAGE \
"===================================================================\n" \
"NLtestMexWrapper: Get nonlinear system for test cases\n" \
"\n" \
"USAGE:\n" \
"  - Constructors:\n" \
"    NLtestMexWrapper( 'new' );\n" \
"n" \
"  - Destructor:\n" \
"    NLtestMexWrapper( 'delete' );\n" \
"n" \
"  - Methods:\n" \
"    NLtestMexWrapper( 'select', name_or_number );\n" \
"    n      = NLtestMexWrapper( 'numberOfTests' );\n" \
"    names  = NLtestMexWrapper( 'listall' );\n" \
"    F      = NLtestMexWrapper( 'evalF', ntest, x [,k] );\n" \
"    JF     = NLtestMexWrapper( 'evalJF', ntest, x );\n" \
"    P      = NLtestMexWrapper( 'pettern', ntest );\n" \
"    n      = NLtestMexWrapper( 'neq', ntest );\n" \
"    ng     = NLtestMexWrapper( 'numGuess', ntest );\n" \
"    g      = NLtestMexWrapper( 'guess', ntest, n );\n" \
"    ne     = NLtestMexWrapper( 'numExact', ntest );\n" \
"    e      = NLtestMexWrapper( 'exact', ntest, n );\n" \
"    ok     = NLtestMexWrapper( 'check', ntest, x );\n" \
"    [L,U]  = NLtestMexWrapper( 'bbox', ntest );\n" \
"    name   = NLtestMexWrapper( 'name', ntest );\n" \
"    bibtex = NLtestMexWrapper( 'bibtex', ntest );\n" \
"\n" \
"===================================================================\n"

using namespace std;

namespace NLproblem {

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_new(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "NLtestMexWrapper('new'): "
    MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    initProblems();
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_delete(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "NLtestMexWrapper('delete'): "
    MEX_ASSERT( nlhs == 0, CMD "expected no output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    // nothing to do
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_numberOfTests(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "NLtestMexWrapper('numberOfTests'): "
    MEX_ASSERT( nlhs == 1, CMD "expected no output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 1, CMD "expected 1 input, nrhs = " << nrhs );
    setScalarInt( arg_out_0, theProblems.size() );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_select(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "NLtestMexWrapper('select',name_or_number): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    integer n;
    if ( mxIsChar(arg_in_1) ) {
      string name = mxArrayToString(arg_in_1);
      try {
        n = theProblemsMap.at(name);
      } catch ( ... ) {
        n = -1;
      }
      ++n;
      MEX_ASSERT(
        n > 0 && n <= integer(theProblems.size()),
        CMD "name = " << name << " cant find"
      );
    } else {
      n = integer( getInt( arg_in_1, CMD "error in get name_or_number" ) );
      MEX_ASSERT(
        n > 0 && n <= integer(theProblems.size()),
        CMD "number = " << n << " out of range"
      );
    }

    mexPrintf( "Selected test N.%d:%s\n", n, theProblems[n-1]->title().c_str() );
    setScalarInt( arg_out_0, n );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_listall(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('listall'): "
    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    integer nprb = integer(theProblems.size());
    arg_out_0 = mxCreateCellMatrix( 1, nprb );

    // Fill cell matrix with input arguments
    for( integer i = 0;  i < nprb; ++i )
      mxSetCell( arg_out_0, i, mxCreateString( theProblems[i]->title().c_str() ) );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_evalF(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('evalF',ntest,...): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];

    #define CMD "NLtestMexWrapper('evalF',ntest,x[,k]): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs >= 2, CMD "expected 3 or 4 input, nrhs = " << nrhs );

    mwSize dimx;
    real_type const * x = getVectorPointer( arg_in_2, dimx, CMD "error in reading x" );

    MEX_ASSERT( dimx == PRB->numEqns(), CMD "bad size(x) = " << dimx );

    dvec_t X(dimx);
    std::copy_n( x, dimx, X.data() );

    if ( nrhs == 3 ) {
      dvec_t F(dimx);
      PRB -> evalF( X, F );
      real_type * f = createMatrixValue( arg_out_0, dimx, 1 );
      std::copy_n( F.data(), dimx, f );
    } else if ( nrhs == 4 ) {
      integer  k  = integer( getInt( arg_in_3, CMD "error in reading k") );
      real_type Fk = PRB -> evalFk( X, k );
      setScalarValue( arg_out_0, Fk );
    } else {
      MEX_ASSERT( false, CMD "expected 3 or 4 input, found " << nrhs );
    }
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_pattern(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('pattern',ntest): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );

    nonlinearSystem * PRB = theProblems[ntest-1];

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    integer n   = PRB->numEqns();
    integer nnz = PRB->jacobianNnz();

    mxArray *args[5];

    real_type * I = createMatrixValue( args[0], 1, nnz );
    real_type * J = createMatrixValue( args[1], 1, nnz );
    real_type * V = createMatrixValue( args[2], 1, nnz );
    setScalarValue( args[3], n );
    setScalarValue( args[4], n );

    ivec_t II(nnz), JJ(nnz);
    PRB->jacobianPattern( II, JJ );

    for ( integer i = 0; i < nnz; ++i ) {
      integer ii = I[i] = II(i)+1; // C to FORTRAN address
      integer jj = J[i] = JJ(i)+1;
      MEX_ASSERT(
        ii > 0 && ii <= n &&  jj > 0 && jj <=  n,
        CMD "idx = " << i << " (i,j) = (" << ii << "," << jj << ") out of range"
      );
      V[i] = 1;
    }

    int ok = mexCallMATLAB( nlhs, plhs, 5, args, "sparse" );

    MEX_ASSERT( ok == 0, CMD "failed the call sparse(...)" );

    #undef CMD

  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_evalJF(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('evalJF',ntest,...): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];

    #define CMD "NLtestMexWrapper('evalJF',ntest,x): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

    mwSize dimx;
    real_type const * x = getVectorPointer( arg_in_2, dimx, CMD "error in reading x" );

    dvec_t X( dimx );
    std::copy_n( x, dimx, X.data() );

    MEX_ASSERT( dimx == PRB->numEqns(), CMD "bad size(x) = " << dimx );

    integer n   = PRB->numEqns();
    integer nnz = PRB->jacobianNnz();

    mxArray *args[5];

    real_type * I = createMatrixValue( args[0], 1, nnz );
    real_type * J = createMatrixValue( args[1], 1, nnz );
    real_type * V = createMatrixValue( args[2], 1, nnz );
    setScalarValue( args[3], n );
    setScalarValue( args[4], n );

    ivec_t II(nnz), JJ(nnz);
    dvec_t VV(nnz);
    PRB->jacobianPattern( II, JJ );
    PRB->jacobian( X, VV );

    for ( integer i = 0; i < nnz; ++i ) {
      integer ii = I[i] = II(i)+1; // C to FORTRAN address
      integer jj = J[i] = JJ(i)+1;
      MEX_ASSERT(
        ii > 0 && ii <= n &&  jj > 0 && jj <=  n,
        CMD "idx = " << i << " (i,j) = (" << ii << "," << jj << ") out of range"
      );
      V[i] = VV(i);
    }

    int ok = mexCallMATLAB( nlhs, plhs, 5, args, "sparse" );

    MEX_ASSERT( ok == 0, CMD "failed the call sparse(...)" );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_neq(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('evalJF',ntest,...): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];

    #define CMD "NLtestMexWrapper('neq',ntest): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    setScalarInt( arg_out_0, PRB->numEqns() );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_numGuess(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('numGuess',ntest): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];
    setScalarInt( arg_out_0, PRB->numInitialPoint() );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_guess(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('guess',ntest,...): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];

    #define CMD "NLtestMexWrapper('guess',ntest,n): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

    integer idx = integer( getInt( arg_in_2, CMD "error in reading n" ) );

    MEX_ASSERT(
      idx > 0 && idx <= PRB->numInitialPoint(),
      CMD "n = " << idx << " out of range"
    );

    dvec_t X( PRB->numEqns() );
    PRB->getInitialPoint( X, idx-1 );
    real_type * x = createMatrixValue( arg_out_0, PRB->numEqns(), 1 );
    std::copy_n( X.data(), PRB->numEqns(), x );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_numExact(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('numExact',ntest): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];
    setScalarInt( arg_out_0, PRB->numExactSolution() );
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_exact(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('exact',ntest,...): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];

    #define CMD "NLtestMexWrapper('exact',ntest,n): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

    integer idx  = integer( getInt( arg_in_2, CMD "error in reading n") );

    MEX_ASSERT(
      idx > 0 && idx <= PRB->numExactSolution(),
      CMD "n = " << idx << " out of range"
    );

    dvec_t X( PRB->numEqns() );
    PRB->getExactSolution( X, idx-1 );
    real_type * x = createMatrixValue( arg_out_0, PRB->numEqns(), 1 );
    std::copy_n( X.data(), PRB->numEqns(), x );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_bbox(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('bbox',ntest): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );

    nonlinearSystem * PRB = theProblems[ntest-1];

    MEX_ASSERT( nlhs == 2, CMD "expected 2 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );

    real_type * L = createMatrixValue( arg_out_0, PRB->numEqns(), 1 );
    real_type * U = createMatrixValue( arg_out_1, PRB->numEqns(), 1 );

    dvec_t LL(PRB->numEqns()), UU(PRB->numEqns());
    PRB->boundingBox( LL, UU );
    std::copy_n( LL.data(), PRB->numEqns(), L );
    std::copy_n( UU.data(), PRB->numEqns(), U );

    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_check(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('check',ntest,...): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );
    #undef CMD

    nonlinearSystem * PRB = theProblems[ntest-1];

    #define CMD "NLtestMexWrapper('check',ntest,x): "

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 3, CMD "expected 3 input, nrhs = " << nrhs );

    mwSize size;
    real_type const * x = getVectorPointer( arg_in_2, size, CMD "error in reading n");

    MEX_ASSERT(
      size == PRB->numEqns(),
      CMD "size(x) = " << size << " expected = " << PRB->numEqns()
    );

    bool ok = true;
    try {
      dvec_t X( PRB->numEqns() );
      std::copy_n( x, PRB->numEqns(), X.data() );
      PRB->checkIfAdmissible( X );
    }
    catch ( runtime_error & err ) {
      //mexPrintf("%s\n",err.what());
      ok = false;
    }
    catch ( ... ) {
      ok = false;
    }
    setScalarBool( arg_out_0, ok );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_name(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    #define CMD "NLtestMexWrapper('name',ntest): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );

    nonlinearSystem * PRB = theProblems[ntest-1];

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input, nrhs = " << nrhs );
    arg_out_0 = mxCreateString( PRB->title().c_str() );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  static
  void
  do_bibtex(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {
    #define CMD "NLtestMexWrapper('bibtex',ntest): "
    integer ntest = integer( getInt( arg_in_1, CMD "error in reading ntest ") );

    MEX_ASSERT(
      ntest > 0 && ntest <= integer(theProblems.size()),
      CMD "ntest = " << ntest << " out of range"
    );

    nonlinearSystem * PRB = theProblems[ntest-1];

    MEX_ASSERT( nlhs == 1, CMD "expected 1 output, nlhs = " << nlhs );
    MEX_ASSERT( nrhs == 2, CMD "expected 2 input"  );
    arg_out_0 = mxCreateString( PRB->bibtex().c_str() );
    #undef CMD
  }

  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

  typedef void (*DO_CMD)( int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] );

  static std::map<std::string,DO_CMD> cmd_to_fun = {
    {"new",do_new},
    {"delete",do_delete},
    {"numberOfTests",do_numberOfTests},
    {"select",do_select},
    {"listall",do_listall},
    {"evalF",do_evalF},
    {"pattern",do_pattern},
    {"evalJF",do_evalJF},
    {"neq",do_neq},
    {"numGuess",do_numGuess},
    {"guess",do_guess},
    {"numExact",do_numExact},
    {"exact",do_exact},
    {"bbox",do_bbox},
    {"check",do_check},
    {"name",do_name},
    {"bibtex",do_bibtex}
  };

  extern "C"
  void
  mexFunction(
    int nlhs, mxArray       *plhs[],
    int nrhs, mxArray const *prhs[]
  ) {

    // the first argument must be a string
    if ( nrhs == 0 ) {
      mexErrMsgTxt(MEX_ERROR_MESSAGE);
      return;
    }

    try {
      MEX_ASSERT(
        mxIsChar(arg_in_0),
        "NLtestMexWrapper: First argument must be a string"
      );
      string cmd = mxArrayToString(arg_in_0);
      DO_CMD pfun = cmd_to_fun.at(cmd);
      pfun( nlhs, plhs, nrhs, prhs );
    } catch ( std::exception const & e ) {
      mexErrMsgTxt( ( string("NLtestMexWrapper Error: ")+e.what() ).c_str() );
    } catch (...) {
  	  mexErrMsgTxt("NLtestMexWrapper failed\n");
    }
  }
}


