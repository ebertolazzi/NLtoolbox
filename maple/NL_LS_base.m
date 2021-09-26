classdef NL_LS_base < handle

  properties (SetAccess = private, Hidden = true)
  end

  methods (Access = protected, Static = true)
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function c = min_quadratic( A, B, do_bound )
      %
      % input stuctures
      % A.t, A.f, A.Df
      % B.t, B.f, B.Df
      %
      %
      % Build cubic P(x) = A*(x-a)^2+B*(x-a)+C such that
      %
      % P(a) = fa, P(b) = fb, P'(a) = Dfa  
      %
      % so that
      %
      % A = (Dfa - (fb - fa)/(b-a))/(a-b)
      % B = Dfa
      % C = fa
      %
      % Search stationary points P'(x) = 2*A*(x-a)+B = 0
      %
      % Setting x-a = r*(b-a)
      %
      % P'(r) = 2*(-Dfa + (fb-fa)/(b-a))*r + Dfa
      %
      % so that the root is
      %
      % r = Dfa/(Dfa - (fb-fa)/(b-a))/2
      %
      r = A.Df/(2*(A.Df - (B.f-A.f)/(B.t-A.t)));
      if do_bound
        r_min = 0.05;
        r_max = 0.95;
        r     = max(r_min,min(r_max,r));
      end
      c = A.t + r*(B.t-A.t);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function c = min_quadratic2( A, B, do_bound )
      r = A.Df/(A.Df - B.Df);
      if do_bound
        r_min = 0.05;
        r_max = 0.95;
        r     = max(r_min,min(r_max,r));
      end
      c = A.t + r*(B.t-A.t);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [c,mono] = min_cubic( A, B, do_bound )
      %
      % Build cubic P(x) = A*(x-a)^3+B*(x-a)^2+C*(x-a)+D such that
      %
      % P(a) = fa, P(b) = fb, P'(a) = Dfa, P'(b) = Dfb  
      %
      % so that
      %
      % A = (Dfa + theta)/(a-b)
      % B = (Dfa + Dfb + 2*theta)/(3*(a-b)^2)
      % C = Dfa
      % D = fa
      %
      theta = (B.Df + A.Df) - 3*(B.f-A.f)/(B.t-A.t);
      %
      % Search stationary points P'(x) = 3*A*(x-a)^2+2*B*(x-a)+C = 0
      %
      % Setting x-a = r*(b-a)
      %
      % P'(r) = (Dfa + Dfb + 2*theta)*r^2 - 2*(Dfa+theta)*r + Dfa
      %
      % the roots are
      %
      % r^+ = ( (Dfa + theta) + sqrt(theta^2-Dfa*Dfb) )/(Dfa + Dfb + 2*theta)
      % r^- = ( -(Dfa + theta) + sqrt(theta^2-Dfa*Dfb) )/(Dfa + Dfb + 2*theta)
      %
      % P''(r^+) = 2*sqrt(theta^2-Dfa*Dfb)/(b-a)
      % P''(r^-) = -2*sqrt(theta^2-Dfa*Dfb)/(b-a)
      %
      % thus the minimum is r^+ for b>a and r^- for b<a
      %
      s     = max(abs([theta,B.Df,A.Df]));
      delta = (theta/s)^2 - (B.Df/s)*(A.Df/s);
      tmp   = A.Df + theta; % -A*(b-a)
      mono  = delta < 0;
      if mono % monotone cubic
        if (B.t-A.t)*A.Df < 0
          r = 1; % minimum at b
        else
          r = 0; % minimum at a
        end
      else
        %
        delta = s*sqrt( delta );
        if B.t>A.t
          % use
          % r^+ = ( (Dfa + theta) + sqrt(theta^2-Dfa*Dfb) )/(Dfa + Dfb + 2*theta)
          if tmp >= 0
            r = (tmp+delta)/(tmp+B.Df+theta);
          else
            % r^+ = Dfa/(Dfa + theta - sqrt(-Dfa*Dfb + theta^2))
            r = A.Df/(tmp-delta);
          end
        else
          % r^- = ( -(Dfa + theta) + sqrt(theta^2-Dfa*Dfb) )/(Dfa + Dfb + 2*theta)
          if tmp <= 0
            r = (delta-tmp)/(tmp+B.Df+theta);
          else
            % r^- = Dfa/(Dfa + theta + sqrt(-Dfa*Dfb + theta^2))
            r = A.Df/(tmp+delta);
          end
        end
      end
      if do_bound
        r_min = 0.05;
        r_max = 0.95;
        r     = max(r_min,min(r_max,r));
      end
      c = A.t + r*(B.t-A.t);
    end
  end
end
