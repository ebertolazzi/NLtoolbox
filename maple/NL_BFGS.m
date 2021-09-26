classdef NL_BFGS < handle

  properties (SetAccess = private, Hidden = true)
    gamma;
    s;
    y;
    rho;
    dim;
    m;
    max_len;
    x0;
    g0;
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = NL_BFGS( varargin ) % dim, nvec
      if nargin == 2
        dim  = varargin{1};
        nvec = varargin{2};
      elseif nargin == 0
        dim  = 0;
        nvec = 0;
      else
        error('NL_BFGS(...) bad number of arguments [%d]\n',nargin);
      end
      self.dim     = dim;
      self.max_len = nvec;
      self.s       = zeros(dim,nvec);
      self.y       = zeros(dim,nvec);
      self.rho     = zeros(1,nvec);
      self.m       = 0;
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
    function setup( self, dim, nvec )
      % Initialize params
      self.dim     = dim;
      self.max_len = nvec;
      self.s       = zeros(dim,nvec);
      self.y       = zeros(dim,nvec);
      self.rho     = zeros(1,nvec);
      self.m       = 0;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function init( self, x0, g0 )
      self.s   = zeros(self.dim,self.max_len);
      self.y   = zeros(self.dim,self.max_len);
      self.rho = zeros(1,self.max_len);
      self.m   = 0;
      self.x0  = x0;
      self.g0  = g0;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function update( self, x1, g1 )
      y   = g1-self.g0;
      s   = x1-self.x0;
      rho = dot(s,y);
      if rho > 0
        yy = dot(y,y);
        self.gamma = rho/yy;
        if self.m >= self.max_len
          % shift
          self.rho(1:end-1) = self.rho(2:end);
          self.s(:,1:end-1) = self.s(:,2:end);
          self.y(:,1:end-1) = self.y(:,2:end);
        else
          self.m = self.m+1;
        end
        self.rho(self.m) = rho;
        self.s(:,self.m) = s(:);
        self.y(:,self.m) = y(:);
      end
      self.x0 = x1;
      self.g0 = g1;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function q = apply( self, x )
      % H(k+1) = (I-s(k)*y(k)^T/rho(k))*H(k)*(I-y(k)*s(k)^T/rho(k))+s(k)*s(k)^T/rho(k)
      q = x;
      if self.m > 0
        alpha = zeros(1,self.m);
        for i=self.m:-1:1
          alpha(i) = dot(self.s(:,i),q)/self.rho(i);
          q        = q - alpha(i)*self.y(:,i);
        end
        % Scale
        q = self.gamma*q;
        for i=1:self.m
          beta = dot(self.y(:,i),q)/self.rho(i);
          q    = q + (alpha(i)-beta)*self.s(:,i);
        end
      end
    end
  end
end
