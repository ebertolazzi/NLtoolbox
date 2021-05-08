classdef NLtest < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    %>
    %> The number of the test used in comput
    %>
    activetest;
  end

  methods
    %
    %> Create a new class instance for the NLtest object.
    %> A class storing \f$ \mathbf{f}(\mathbf{x}) \f$ and its jacobian, guess, 
    %> and exact solution (if available) for many nonlinear systems.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    ref = NLtest()
    %>
    %> \endrst
    %>
    %> **Output**
    %>
    %> - `ref` the instance of the Matlab class.
    %>
    function self = NLtest()
      NLtestMexWrapper( 'new' );
      self.activetest = 0 ;
    end

    function delete(self)
      NLtestMexWrapper( 'delete' );
    end

    %>
    %> Print on console the list of available test cases.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    ref.listall();
    %>
    %> \endrst
    %>
    function names = listall( self )
      names = NLtestMexWrapper( 'listall' );
    end

    %>
    %> Return the number of available test cases.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    ntest = ref.numTest();
    %>
    %> \endrst
    %>
    function n = numTest( self )
      n = NLtestMexWrapper( 'numTest' );
    end

    %>
    %> Select the active test to be used.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    ntest = 23;
    %>    ref.select( ntest );
    %>
    %>    test_name = 'Beale';
    %>    ref.select( test_name );
    %>
    %> \endrst
    %>
    function select( self, name_or_number )
      self.activetest = NLtestMexWrapper( 'select', name_or_number );
    end

    %>
    %> Evaluate \f$ \mathbf{f}(\mathbf{x}) \f$ or \f$ f_k(\mathbf{x}) \f$ for the active test.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    f = ref.evalF( x );
    %>
    %>    k = 2;
    %>    f2 = ref.select( x, k ); % Only k-th component
    %>
    %> \endrst
    %>
    function F = evalF( self, varargin )
      % varargin = {x} 
      % varargin = {x, ncomp}
      F = NLtestMexWrapper( 'evalF', self.activetest, varargin{:} );
    end

    %>
    %> Return a sparse matrix containing
    %>
    %> \f[ \dfrac{\partial\mathbf{f}(\mathbf{x})}{\partial\mathbf{x}} \f]
    %>
    %> for the active test.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    J = ref.evalJF( x );
    %>
    %> \endrst
    %>
    function JF = evalJF( self, x )
      JF = NLtestMexWrapper( 'evalJF', self.activetest, x );
    end

    %>
    %> Return a sparse matrix containing the pattern of nonzeros of the jacobian
    %>
    %> \f[ \dfrac{\partial\mathbf{f}(\mathbf{x})}{\partial\mathbf{x}} \f]
    %>
    %> for the active test.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    Jpattern = ref.pattern( x );
    %>
    %> \endrst
    %>
    function P = pattern( self )
      P = NLtestMexWrapper( 'pattern', self.activetest );
    end

    %>
    %> Return number of equations of the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    N = ref.neq();
    %>
    %> \endrst
    %>
    function n = neq( self )
      n = NLtestMexWrapper( 'neq', self.activetest );
    end

    %>
    %> Return number of initial guess points of the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    NG = ref.numGuess();
    %>
    %> \endrst
    %>
    function n = numGuess( self )
      n = NLtestMexWrapper( 'numGuess', self.activetest );
    end

    %>
    %> Return the n-th guess point of the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    g = ref.guess( n );
    %>
    %> \endrst
    %>
    function g = guess( self, n )
      g = NLtestMexWrapper( 'guess', self.activetest, n );
    end

    %>
    %> Return number of available exact solution of the active nonlinear system.
    %> May be 0.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    NG = ref.numExact();
    %>
    %> \endrst
    %>
    function n = numExact( self )
      n = NLtestMexWrapper( 'numExact', self.activetest );
    end

    %>
    %> Return the n-th exact solution of the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    e = ref.exact( n );
    %>
    %> \endrst
    %>
    function sol = exact( self, n )
      sol = NLtestMexWrapper( 'exact', self.activetest, n );
    end

    %>
    %> Check if `x` satisfy the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    ok = ref.check( x );
    %>
    %> \endrst
    %>
    function ok = check( self, x )
      ok = NLtestMexWrapper( 'check', self.activetest, x );
    end

    %>
    %> Return lower and upper bound for the solutions of the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    [L,U] = ref.bbox( x );
    %>
    %> \endrst
    %>
    function [L,U] = bbox( self )
      [L,U] = NLtestMexWrapper( 'bbox', self.activetest );
    end

    %>
    %> Return the number of the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    n = ref.numActive( x );
    %>
    %> \endrst
    %>
    function n = numActive( self )
      n = self.activetest;
    end

    %>
    %> Return the name (string) of the n-th or the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    testname = ref.name( 10 ); % the name of the test n. 10
    %>
    %>    testname = ref.name(); % the name of the active test
    %>
    %> \endrst
    %>
    function testname = name( self, varargin )
      if nargin > 1
        testname = NLtestMexWrapper( 'name', varargin{1} );
      else
        testname = NLtestMexWrapper( 'name', self.activetest );
      end
    end

    %>
    %> Return the bibtex string of the n-th or the active nonlinear system.
    %>
    %> **Usage:**
    %>
    %> \rst
    %>
    %> .. code-block:: matlab
    %>
    %>    bib = ref.bibtex( 10 ); % the bibtex of the test n. 10
    %>
    %>    bib = ref.bibtex(); % the bibtex of the active test
    %>
    %> \endrst
    %>
    function bib = bibtex( self, varargin )
      if nargin > 1
        bib = NLtestMexWrapper( 'bibtex', varargin{1} );
      else
        bib = NLtestMexWrapper( 'bibtex', self.activetest );
      end
    end

  end
end