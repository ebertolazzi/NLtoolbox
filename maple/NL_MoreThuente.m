%CVSRCH   More-Thuente line search from MINPACK.
%
% Re-factor of API for use in the Poblano Toolbox, 
% Daniel M. Dunlavy, August 2008, March 2009
%    
% Translation of minpack subroutine cvsrch
% Dianne O'Leary   July 1991
% **********
%
% Subroutine cvsrch
%
% The purpose of cvsrch is to find a step which satisfies 
% a sufficient decrease condition and a curvature condition.
% The user must provide a subroutine which calculates the
% function and the gradient.
%
% At each stage the subroutine updates an interval of
% uncertainty with endpoints stx and sty. The interval of
% uncertainty is initially chosen so that it contains a 
% minimizer of the modified function
%
%      f(x+stp*s) - f(x) - f_tolerance*stp*(gradf(x)'s).
%
% If a step is obtained for which the modified function 
% has a nonpositive function value and nonnegative derivative, 
% then the interval of uncertainty is chosen so that it 
% contains a minimizer of f(x+stp*s).
%
% The algorithm is designed to find a step which satisfies 
% the sufficient decrease condition 
%
%           f(x+stp*s) <= f(x) + f_tolerance*stp*(gradf(x)'s),
%
% and the curvature condition
%
%           abs(gradf(x+stp*s)'s)) <= g_tolerance*abs(gradf(x)'s).
%
% If f_tolerance is less than g_tolerance and if, for example, the function
% is bounded below, then there is always a step which satisfies
% both conditions. If no step can be found which satisfies both
% conditions, then the algorithm usually stops when rounding
% errors prevent further progress. In this case stp only 
% satisfies the sufficient decrease condition.
%
% The subroutine statement is
%
%        subroutine cvsrch(fcn,n,x,f,g,s,stp,f_tolerance,g_tolerance,x_tolerance,
%                          step_min,step_max,maxfev,info,nfev,wa)
% where
%
% fcn is the name of the user-supplied subroutine which
% calculates the function and the gradient.  fcn must 
% be declared in an external statement in the user 
% calling program, and should be written as follows.
%
%	  subroutine fcn(n,x,f,g)
%         integer n
%         f
%         x(n),g(n)
%	  ----------
%         Calculate the function at x and
%         return this value in the variable f.
%         Calculate the gradient at x and
%         return this vector in g.
%	  ----------
%	  return
%	  end
%
%       n is a positive integer input variable set to the number
%	  of variables.
%
%	x is an array of length n. On input it must contain the
%	  base point for the line search. On output it contains 
%         x + stp*s.
%
%	f is a variable. On input it must contain the value of f
%         at x. On output it contains the value of f at x + stp*s.
%
%	g is an array of length n. On input it must contain the
%         gradient of f at x. On output it contains the gradient
%         of f at x + stp*s.
%
%	s is an input array of length n which specifies the
%         search direction.
%
%	stp is a nonnegative variable. On input stp contains an
%         initial estimate of a satisfactory step. On output
%         stp contains the final estimate.
%
%       f_tolerance and g_tolerance are nonnegative input variables. Termination
%         occurs when the sufficient decrease condition and the
%         directional derivative condition are satisfied.
%
%	x_tolerance is a nonnegative input variable. Termination occurs
%         when the relative width of the interval of uncertainty 
%	  is at most x_tolerance.
%
%	step_min and step_max are nonnegative input variables which 
%	  specify lower and upper bounds for the step.
%
%	maxfev is a positive integer input variable. Termination
%         occurs when the number of calls to fcn is at least
%         maxfev by the end of an iteration.
%
%	info is an integer output variable set as follows:
%
%	  info = 0  Improper input parameters.
%
%	  info = 1  The sufficient decrease condition and the
%                   directional derivative condition hold.
%
%	  info = 2  Relative width of the interval of uncertainty
%		    is at most x_tolerance.
%
%	  info = 3  Number of calls to fcn has reached maxfev.
%
%	  info = 4  The step is at the lower bound step_min.
%
%	  info = 5  The step is at the upper bound step_max.
%
%	  info = 6  Rounding errors prevent further progress.
%                   There may not be a step which satisfies the
%                   sufficient decrease and curvature conditions.
%                   Tolerances may be too small.
%
%       nfev is an integer output variable set to the number of
%         calls to fcn.
%
%	wa is a work array of length n.
%
% Argonne National Laboratory. MINPACK Project. June 1983
% Jorge J. More', David J. Thuente
%
% **********


%CSTEP   More-Thuente line search step from MINPACK.
%
% Translation of minpack subroutine cstep 
% Dianne O'Leary   July 1991
% **********
%
% Subroutine cstep
%
% The purpose of cstep is to compute a safeguarded step for
% a linesearch and to update an interval of uncertainty for
% a minimizer of the function.
%
% The parameter stepx.t contains the step with the least function value.
% The parameter stepp.t contains the current step. It is assumed that the
% derivative at stepx.t is negative in the direction of the step.
% If brackt is set true then a minimizer has been bracketed in an interval
% of uncertainty with endpoints stepx.t and stepy.t.
%
% The subroutine statement is
%
%   cstep(stepx.[t,f,Df],stepy.[t,f,Df],stepp.[t,f,Df],brackt,step_min,step_max,info)
%
% where
%
% stepx.[t,f,Df] are variables which specify the step,
%   the function, and the derivative at the best step obtained
%   so far. The derivative must be negative in the direction
%   of the step, that is, stepx.Df and stepp.t-stepx.t must have opposite 
%   signs. On output these parameters are updated appropriately.
%
% stepy.[t,f,Df] are variables which specify the step,
%   the function, and the derivative at the other endpoint of
%   the interval of uncertainty. On output these parameters are 
%   updated appropriately.
%
% stepp.[t,f,Df] are variables which specify the step,
%   the function, and the derivative at the current step.
%   If brackt is set true then on input stepp.t must be
%   between stepx.t and stepy.t. On output stepp.t is set to the new step.
%
% brackt is a logical variable which specifies if a minimizer
%   has been bracketed. If the minimizer has not been bracketed
%   then on input brackt must be set false. If the minimizer
%   is bracketed then on output brackt is set true.
%
% step_min and step_max are input variables which specify lower 
%   and upper bounds for the step.
%
% info is an integer output variable set as follows:
%   If info = 1,2,3,4,5, then the step has been computed
%   according to one of the five cases below. Otherwise
%   info = 0, and this indicates improper input parameters.
%
% Argonne National Laboratory. MINPACK Project. June 1983
% Jorge J. More', David J. Thuente
%
% **********
% Redistribution and use in source and binary forms, with or
% without modification, are permitted provided that the
% following conditions are met:
% 
% 1. Redistributions of source code must retain the above
% copyright notice, this list of conditions and the following
% disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above
% copyright notice, this list of conditions and the following
% disclaimer in the documentation and/or other materials
% provided with the distribution.
% 
% 3. The end-user documentation included with the
% redistribution, if any, must include the following
% acknowledgment:
% 
%    "This product includes software developed by the
%    University of Chicago, as Operator of Argonne National
%    Laboratory.
% 
% Alternately, this acknowledgment may appear in the software
% itself, if and wherever such third-party acknowledgments
% normally appear.
% 
% 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
% WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
% UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
% THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
% OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
% OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
% USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
% THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
% DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
% UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
% BE CORRECTED.
% 
% 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
% HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
% ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
% INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
% ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
% PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
% SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
% (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
% EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
% POSSIBILITY OF SUCH LOSS OR DAMAGES.
% 
% Redistribution and use in source and binary forms, with or
% without modification, are permitted provided that the
% following conditions are met:
% 
% 1. Redistributions of source code must retain the above
% copyright notice, this list of conditions and the following
% disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above
% copyright notice, this list of conditions and the following
% disclaimer in the documentation and/or other materials
% provided with the distribution.
% 
% 3. The end-user documentation included with the
% redistribution, if any, must include the following
% acknowledgment:
% 
%    "This product includes software developed by the
%    University of Chicago, as Operator of Argonne National
%    Laboratory.
% 
% Alternately, this acknowledgment may appear in the software
% itself, if and wherever such third-party acknowledgments
% normally appear.
% 
% 4. WARRANTY DISCLAIMER. THE SOFTWARE IS SUPPLIED "AS IS"
% WITHOUT WARRANTY OF ANY KIND. THE COPYRIGHT HOLDER, THE
% UNITED STATES, THE UNITED STATES DEPARTMENT OF ENERGY, AND
% THEIR EMPLOYEES: (1) DISCLAIM ANY WARRANTIES, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
% OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE
% OR NON-INFRINGEMENT, (2) DO NOT ASSUME ANY LEGAL LIABILITY
% OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR
% USEFULNESS OF THE SOFTWARE, (3) DO NOT REPRESENT THAT USE OF
% THE SOFTWARE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS, (4)
% DO NOT WARRANT THAT THE SOFTWARE WILL FUNCTION
% UNINTERRUPTED, THAT IT IS ERROR-FREE OR THAT ANY ERRORS WILL
% BE CORRECTED.
% 
% 5. LIMITATION OF LIABILITY. IN NO EVENT WILL THE COPYRIGHT
% HOLDER, THE UNITED STATES, THE UNITED STATES DEPARTMENT OF
% ENERGY, OR THEIR EMPLOYEES: BE LIABLE FOR ANY INDIRECT,
% INCIDENTAL, CONSEQUENTIAL, SPECIAL OR PUNITIVE DAMAGES OF
% ANY KIND OR NATURE, INCLUDING BUT NOT LIMITED TO LOSS OF
% PROFITS OR LOSS OF DATA, FOR ANY REASON WHATSOEVER, WHETHER
% SUCH LIABILITY IS ASSERTED ON THE BASIS OF CONTRACT, TORT
% (INCLUDING NEGLIGENCE OR STRICT LIABILITY), OR OTHERWISE,
% EVEN IF ANY OF SAID PARTIES HAS BEEN WARNED OF THE
% POSSIBILITY OF SUCH LOSS OR DAMAGES.
% 

classdef NL_MoreThuente < NL_LS_base

  properties (SetAccess = private, Hidden = true)
    x_tolerance;
    f_tolerance;
    g_tolerance;
    step_min;
    step_max;
    max_iter;
    iter;
  end

  methods (Access = public, Static = true, Hidden = false)
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function str = into_to_string( info )
      switch info
      case 0; str = 'Improper input parameters.';
      case 1; str = [ 'The sufficient decrease condition and the\n', ...
                      'directional derivative condition hold.'];
      case 2; str = 'Relative width of the interval of uncertainty is at most x_tolerance.';
      case 3; str = 'Number of calls to fcn has reached maxfev.';
      case 4; str = 'The step is at the lower bound step_min.';
      case 5; str = 'The step is at the upper bound step_max.';
      case 6; str = [ 'Rounding errors prevent further progress.\n', ...
                      'There may not be a step which satisfies the\n', ...
                      'sufficient decrease and curvature conditions.\n', ...
                      'Tolerances may be too small.'];
      otherwise
        str = 'Unknown error';
      end
    end
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = NL_MoreThuente( )
      % Initialize params
      self.x_tolerance = 1e-8;
      self.f_tolerance = [0.1,0.01,0.0001];%1e-8;
      self.g_tolerance = [0.1,0.5,0.8];%1e-8;
      self.step_min    = 1e-10;
      self.step_max    = 1e+10;
      self.max_iter    = [10,10,20];
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function delete( ~ )
      %% Destroy the C++ class instance
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function set_x_tolerance( self, x_tolerance )
      if x_tolerance < 0
        error('NL_MoreThuente: set_x_tolerance(%g) argument must be > 0\n',x_tolerance);
      end
      self.x_tolerance = x_tolerance;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function set_f_tolerance( self, f_tolerance )
      if f_tolerance < 0
        error('NL_MoreThuente: set_f_tolerance(%g) argument must be > 0\n',f_tolerance);
      end
      self.f_tolerance = f_tolerance;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function set_g_tolerance( self, g_tolerance )
      if g_tolerance < 0
        error('NL_MoreThuente: set_g_tolerance(%g) argument must be > 0\n',g_tolerance);
      end
      self.g_tolerance = g_tolerance;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function set_minmax( self, step_min, step_max )
      if step_min < 0 || step_max < step_min
        error( ...
          'NL_MoreThuente: set_minmax(min=%g,max=%g) argument must be 0 < min < max\n', ...
          step_min, step_max ...
        );
      end
      self.step_min = step_min;
      self.step_max = step_max;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function iter = get_iter( self )
      iter = self.iter;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [stepx,stepy,stepp,brackt,info] = refine(self,stepx,stepy,stepp,brackt,step_min,step_max)
      info = 0;
      %
      % Check the input parameters for errors.
      %
      mi = min(stepx.t,stepy.t);
      ma = max(stepx.t,stepy.t);
      if brackt && (stepp.t <= mi || stepp.t >= ma)
        return
      end
      if stepx.Df*(stepp.t-stepx.t) >= 0 || step_max < step_min
        return
      end
      %
      % Determine if the derivatives have opposite sign.
      %
      sgnd = stepp.Df*sign(stepx.Df);
      if stepp.f > stepx.f 
        %
        % First case.
        % A higher function value. The minimum is bracketed.
        % If the cubic step is closer to stepx.t than the quadratic step,
        % the cubic step is taken, else the average of the cubic and
        % quadratic steps is taken.
        %
        info   = 1;
        bound  = true;
        brackt = true;
        [stpc,~] = self.min_cubic( stepx, stepp, true );
        stpq     = self.min_quadratic( stepx, stepp, true );
        if abs(stpc-stepx.t) < abs(stpq-stepx.t) 
          stpf = stpc;
        else
          stpf = stpc + (stpq - stpc)/2;
        end
      elseif sgnd < 0
        %
        % Second case.
        % A lower function value and derivatives of opposite sign.
        % The minimum is bracketed.
        % If the cubic step is closer to stepx.t than the quadratic (secant) step, 
        % the cubic step is taken, else the quadratic step is taken.
        %
        info   = 2;
        bound  = false;
        brackt = true;
        [stpc,~] = self.min_cubic( stepp, stepx, true );
        stpq     = self.min_quadratic2( stepp, stepx, true );
        if abs(stpc-stepp.t) > abs(stpq-stepp.t)
          stpf = stpc;
        else
          stpf = stpq;
        end
      elseif abs(stepp.Df) < abs(stepx.Df)
        %
        % Third case.
        % A lower function value, derivatives of the same sign, 
        % and the magnitude of the derivative decreases.
        % The cubic step is only used if
        %   - the cubic tends to infinity in the direction of the step or
        %   - if the minimum of the cubic is beyond stepp.t.
        % Otherwise the cubic step is defined to be either step_min or step_max.
        % The quadratic (secant) step is also computed and if the minimum
        % is bracketed then the the step closest to stepx.t is taken,
        % else the step farthest away is taken.
        %
        info  = 3;
        bound = true;
        [stpc,mono] = self.min_cubic( stepp, stepx, false ); % @@@@@@@@@@@@
        if sign(stepx.t - stepp.t) == sign(stpc - stepp.t) || mono
          if stepp.t > stepx.t
            stpc = self.step_max;
          else
            stpc = self.step_min;
          end
        end
    
        stpq = self.min_quadratic2( stepp, stepx, false );
        if brackt 
          if abs(stepp.t-stpc) < abs(stepp.t-stpq)
            stpf = stpc;
          else
            stpf = stpq;
          end
        else
          if abs(stepp.t-stpc) > abs(stepp.t-stpq)
            stpf = stpc;
          else
            stpf = stpq;
          end 
        end
      else
        %
        % Fourth case. A lower function value, derivatives of the
        % same sign, and the magnitude of the derivative does
        % not decrease. If the minimum is not bracketed, the step
        % is either step_min or step_max, else the cubic step is taken.
        %
        info  = 4;
        bound = false;
        if brackt
          [stpf,~] = self.min_cubic( stepp, stepy, true );
        elseif stepp.t > stepx.t
          stpf = self.step_max;
        else
          stpf = self.step_min;
        end 
      end 
      %
      % Update the interval of uncertainty. This update does not
      % depend on the new step or the case analysis above.
      %
      if stepp.f > stepx.f
        stepy = stepp;
      else
        if sgnd < 0
          stepy = stepx;
        end 
        stepx = stepp;
      end
      %
      % Compute the new step and safeguard it.
      %
      stepp.t = min(step_max,max(step_min,stpf));
      if brackt && bound
        newstp = stepx.t+(2/3)*(stepy.t-stepx.t);
        if stepy.t > stepx.t
          stepp.t = min(newstp,stepp.t);
        else
          stepp.t = max(newstp,stepp.t);
        end
      end
      return
      %
      % Last card of subroutine refine.
      %
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [step,info,nfev] = linesearch(self,fcn,stp,f0,Df0)
      % Initialize params
      
      nfev   = 0;
      xtrapf = 4;
      info   = 0;
      infoc  = 1;
      %
      % Check the input parameters for errors.
      %
      if stp <= 0 || Df0 >= 0
        step.t  = 0;
        step.f  = f0;
        step.Df = Df0;
        return;
      end
      %
      % Initialize local variables.
      %
      brackt = false;
      stage1 = true;
      width  = self.step_max - self.step_min;
      width1 = 2*width;
      %
      % The variables stx, fx, dx contain the values of the step, 
      % function, and directional derivative at the best step.
      % The variables sty, fy, dy contain the value of the step,
      % function, and derivative at the other endpoint of
      % the interval of uncertainty.
      % The variables stp, f, dg contain the values of the step,
      % function, and derivative at the current step.
      %
      stepx.t  = 0;
      stepx.f  = f0;
      stepx.Df = Df0;
      stepy.t  = 0.0;
      stepy.f  = f0;
      stepy.Df = Df0;
      %
      % Start of iteration.
      %
      step.t    = stp;
      self.iter = 0;
      for jjj=1:3
        
        minus_gtol_Df0    = -self.g_tolerance(jjj)*Df0;
        min_ftol_gtol_Df0 = min(self.f_tolerance(jjj),self.g_tolerance(jjj))*Df0;
        dgtest            = self.f_tolerance(jjj)*Df0;

        for kkk=1:self.max_iter(jjj)
          self.iter = self.iter+1;
          %
          % Set the minimum and maximum steps to correspond
          % to the present interval of uncertainty.
          %
          if brackt
            stmin = min(stepx.t,stepy.t);
            stmax = max(stepx.t,stepy.t);
          else
            stmin = stepx.t;
            stmax = step.t + xtrapf*(step.t - stepx.t);
          end 
          %
          % Force the step to be within the bounds step_max and step_min.
          %
          step.t = max(min(step.t,self.step_max),self.step_min);
          %
          % If an unusual termination is to occur then let 
          % stp be the lowest point obtained so far.
          %
          if (brackt && (step.t <= stmin || step.t >= stmax)) || infoc == 0 || ...
              (brackt && stmax-stmin <= self.x_tolerance*stmax)
            step = stepx;
          end
          %
          % Evaluate the function and gradient at step.t
          % and compute the directional derivative.
          %
          [f,Df] = feval(fcn,step.t);
          nfev   = nfev + 1;
          ftest1 = f0 + step.t*dgtest;
    
          step.f  = f;
          step.Df = Df;
          if isnan(f) || isnan(Df)
            fprintf('\n\n\nFound NAN\n\n'); 
            step.f  = Inf;
            step.Df = Inf;           
          end
          %
          % Test for convergence.
          %
          if (brackt && (step.t <= stmin || step.t >= stmax)) || infoc == 0
            info = 6;
          end
          if step.t >= self.step_max && step.f <= ftest1 && step.Df <= dgtest
            info = 5;
          end
          if step.t <= self.step_min && (step.f > ftest1 || step.Df >= dgtest)
            info = 4;
          end
          %if nfev >= self.maxfev 
          %  info = 3;
          %end
          if brackt && stmax-stmin <= self.x_tolerance*stmax
            info = 2;
          end
          if step.f <= ftest1 && abs(step.Df) <= minus_gtol_Df0
            info = 1;
          end
          %
          % Check for termination.
          %
          if info ~= 0
            return;
          end
          %
          % In the first stage we seek a step for which the modified
          % function has a nonpositive value and nonnegative derivative.
          %
          if stage1 && step.f <= ftest1 && step.Df >= min_ftol_gtol_Df0
            stage1 = false;
          end
          %
          % A modified function is used to predict the step only if
          % we have not obtained a step for which the modified
          % function has a nonpositive function value and nonnegative 
          % derivative, and if a lower function value has been  
          % obtained but the decrease is not sufficient.
          %
          if stage1 && step.f <= stepx.f && step.f > ftest1
            %
            % Define the modified function and derivative values.
            %
            step.f  = step.f  - step.t*dgtest;
            step.Df = step.Df - dgtest;

            stepx.f  = stepx.f  - stepx.t*dgtest;
            stepx.Df = stepx.Df - dgtest;

            stepy.f  = stepy.f  - stepy.t*dgtest;
            stepy.Df = stepy.Df - dgtest;
            % 
            % Call refine to update the interval of uncertainty and to compute the new step.
            %
            [stepx,stepy,step,brackt,infoc] = self.refine(stepx,stepy,step,brackt,stmin,stmax);
            %
            % Reset the function and gradient values for f.
            %
            step.f  = step.f  + step.t*dgtest;
            step.Df = step.Df + dgtest;

            stepx.f  = stepx.f  + stepx.t*dgtest;
            stepx.Df = stepx.Df + dgtest;

            stepy.f  = stepy.f  + stepy.t*dgtest;
            stepy.Df = stepy.Df + dgtest;
          else
            %
            % Call refine to update the interval of uncertainty 
            % and to compute the new step.
            %
            [stepx,stepy,step,brackt,infoc] = self.refine(stepx,stepy,step,brackt,stmin,stmax);
          end
          %
          % Force a sufficient decrease in the size of the
          % interval of uncertainty.
          %
          if (brackt) 
            if (abs(stepy.t-stepx.t) >= 0.66*width1) 
              step.t = stepx.t + (stepy.t - stepx.t)/2;
            end
            width1 = width;
            width  = abs(stepy.t-stepx.t);
          end
          %if kkk > 100
          %  fprintf('failure\n'); 
          %end
        end
        %
        % End of iteration.
        %
      end
      info = 3;
      %
      % Last card of subroutine cvsrch.
      %
    end
  end
end
