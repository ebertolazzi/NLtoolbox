classdef NL_CG < handle

  properties (SetAccess = private, Hidden = true)
    old_x;
    old_g;
    old_d;
    do_restart;
    method;
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = NL_CG()
      self.do_restart = true;
      self.method     = 'FR';
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
      self.do_restart = true;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function d = update( self, x, g )
      x0 = self.old_x;
      g0 = self.old_g;
      d0 = self.old_d;
      if ~self.do_restart
        s  = x-x0;
        y  = g-g0;
        switch self.method
        case 'HS' % Hestenes-Stiefel
          bt = dot(g,y)/dot(y,d0);
        case 'FR' % Fletcher-Reeves
          bt = dot(g,g)/dot(g0,g0);
        case 'PRP' % Polak-Ribiere
          bt = dot(g,y)/dot(g0,g0);
        case 'GN' % Gilbert-Nocedal
          g0g0    = dot(g0,g0);
          beta_FR = dot(g,g)/g0g0;
          beta_PR = dot(g,y)/g0g0;
          bt = max(-beta_FR,min(beta_PR,beta_FR));
        case 'CD'
          bt = -dot(g,g)/dot(g0,d0);
        case 'LS'
          bt = -dot(g,y)/dot(g0,d0);
        case 'DY'
          bt = dot(g,g)/dot(y,d0);
        case 'WYL'
          t1 = dot(g,g);
          t2 = dot(d0,d0);
          bt = dot(g,g-sqrt(t1/t2)*g0)/t2;
        case 'RMIL'
          bt = dot(g,y)/dot(d0,d0);
        case 'NPPRP'
          t1 = dot(g,g);
          t2 = dot(g0,g0);
          t3 = dot(d0,d0);
          bt = (t1-sqrt(t1/t3)*abs(dot(g,d0)))/t2;
        case 'WAM'
          t1 = dot(g,g);
          t2 = dot(g0,g0);
          t3 = dot(d0,d0);
          bt = (t1-sqrt(t1/t3)*abs(dot(g0,d0)))/t2;
        otherwise
          error('NLCG: update, unknown method %s',self.method);
        end
        d = -g + bt*d0;
      else
        d = -g;
      end

      self.old_x = x;
      self.old_g = g;
      self.old_d = d;
      self.do_restart = false;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
end
