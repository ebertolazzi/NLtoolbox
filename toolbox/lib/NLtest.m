classdef NLtest < handle
  %% MATLAB class wrapper for the underlying C++ class
  properties (SetAccess = private, Hidden = true)
    activetest;
  end

  methods
    function self = NLtest()
      %% Create a new C++ class instance for the NLtest object
      %
      %  Usage:
      %    ref = NLtest()
      %
      %  On output:
      %    ref: reference handle to the object instance
      %
      NLtestMexWrapper( 'new' );
      self.activetest = 0 ;
    end

    function delete(self)
      %% Destroy the C++ class instance
      NLtestMexWrapper( 'delete' );
    end

    function names = listall( self )
      names = NLtestMexWrapper( 'listall' );
    end

    function n = numTest( self )
      n = NLtestMexWrapper( 'numTest' );
    end

    function select( self, name_or_number )
      self.activetest = NLtestMexWrapper( 'select', name_or_number );
    end

    function F = evalF( self, varargin )
      % varargin = {x} 
      % varargin = {x, ncomp}
      F = NLtestMexWrapper( 'evalF', self.activetest, varargin{:} );
    end

    function JF = evalJF( self, x )
      JF = NLtestMexWrapper( 'evalJF', self.activetest, x );
    end

    function P = pattern( self )
      P = NLtestMexWrapper( 'pattern', self.activetest );
    end

    function n = neq( self )
      n = NLtestMexWrapper( 'neq', self.activetest );
    end

    function n = numGuess( self )
      n = NLtestMexWrapper( 'numGuess', self.activetest );
    end

    function g = guess( self, n )
      g = NLtestMexWrapper( 'guess', self.activetest, n );
    end

    function n = numExact( self )
      n = NLtestMexWrapper( 'numExact', self.activetest );
    end

    function sol = exact( self, n )
      sol = NLtestMexWrapper( 'exact', self.activetest, n );
    end

    function ok = check( self, x )
      ok = NLtestMexWrapper( 'check', self.activetest, x );
    end

    function [L,U] = bbox( self )
      [L,U] = NLtestMexWrapper( 'bbox', self.activetest );
    end

    function n = numActive( self )
      n = self.activetest;
    end

    function testname = name( self )
      testname = NLtestMexWrapper( 'name', self.activetest );
    end

    function testname = info( self )
      testname = NLtestMexWrapper( 'info', self.activetest );
    end

  end
end