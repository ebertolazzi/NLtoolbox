classdef NL_BB < handle

  properties (SetAccess = private, Hidden = true)
    old_f;
    old_x;
    old_g;
    method;
    k_step;
    alpha_stored;
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = NL_BB()
      self.method       = 1;
      self.k_step       = 0;
      self.alpha_stored = ones(1,5);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function delete( ~ )
      %% Destroy the C++ class instance
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function setup( self, method )
      % Initialize params
      self.method = method;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function restart( self )
      self.k_step       = 0;
      self.alpha_stored = ones(1,5);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function alpha = update( self, f, x, g )
      f0 = self.old_f;
      x0 = self.old_x;
      g0 = self.old_g;
      if self.k_step > 0
        s  = x-x0;
        y  = g-g0;
        sy = dot(s,y);
        ss = dot(s,s);
        switch self.method
        case 1
          yy    = dot(y,y);
          alpha = sy/(yy);
          if alpha <= 1e-10 || alpha > 1e10
            alpha = 1;
          end
        case 2
          alpha = ss/sy;
          if alpha <= 1e-10 || alpha > 1e10
            alpha = 1;
          end
        case 3
          % Conic Interpolation ('Modified BB')
          alpha = ss/sy;
          if alpha <= 1e-10 || alpha > 1e10
            alpha = 1;
          end
          alphaConic = ss/(6*(f0 - f) + 4*dot(g,s) + 2*dot(g0,s));
          if alphaConic > 0.001*alpha && alphaConic < 1000*alpha
            alpha = alphaConic;
          end
        case 4
          % Gradient Method with retards (bb type 1, random selection of previous step)
          alpha = ss/sy;
          if alpha <= 1e-10 || alpha > 1e10
            alpha = 1;
          end
          self.alpha_stored(1+mod(self.k_step-2,5)) = alpha;
          alpha = self.alpha_stored(ceil(rand*5));
        otherwise
          error('NL_BB: update, unknown method %s',self.method);
        end
      else
        alpha = 0;
      end

      self.old_f  = f;
      self.old_x  = x;
      self.old_g  = g;
      self.k_step = self.k_step+1;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
end
