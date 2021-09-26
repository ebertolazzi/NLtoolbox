classdef NL_BROYDEN < handle

  properties (SetAccess = private, Hidden = true)
    x0;
    F0;
    B0;
    v;
    w;
    z;
    dim;
    m;
    max_len;
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = NL_BROYDEN( varargin ) % dim, nvec
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
      self.v       = zeros(dim,nvec);
      self.w       = zeros(dim,nvec);
      self.z       = zeros(dim,nvec);
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
      self.v       = zeros(dim,nvec);
      self.w       = zeros(dim,nvec);
      self.z       = zeros(dim,nvec);
      self.m       = 0;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function init( self, x0, F0, B0 )
      self.m  = 0;
      self.x0 = x0;
      self.F0 = F0;
      self.B0 = B0;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function is_f = is_full( self )
      is_f = self.m >= self.max_len;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function update( self, lambda, x1, F1 )
      s      = x1-self.x0;
      slen   = sqrt(dot(s,s));
      z      = F1-(1-lambda)*self.F0;
      zz     = self.apply( z );
      alpha  = slen+dot(zz,s)/slen;
      self.m = self.m+1;
      if self.m > self.max_len
        error('NL_BROYDEN::update, exceeded max length [$d]', self.max_len);
      end
      self.z(:,self.m) = z/slen;
      self.v(:,self.m) = s/slen;
      self.w(:,self.m) = zz/alpha;
      self.x0          = x1;
      self.F0          = F1;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function q = mult( self, x )
      q = self.B0*x;
      if self.m > 0
        for i=1:self.m
          t = dot(x,self.v(:,i));
          q = q + t*self.z(:,i);
        end
      end
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function q = mult_transposed( self, x )
      q = self.B0.'*x;
      if self.m > 0
        for i=1:self.m
          t = dot(x,self.z(:,i));
          q = q + t*self.v(:,i);
        end
      end
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function q = apply( self, x )
      q = lsqr(self.B0,x);
      if self.m > 0
        for i=1:self.m
          t = dot(q,self.v(:,i));
          q = q - t*self.w(:,i);
        end
      end
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function q = apply_transposed( self, x )
      q = x;
      if self.m > 0
        for i=self.m:-1:1
          t = dot(q,self.w(:,i));
          q = q - t*self.v(:,i);
        end
      end
      q = lsqr(self.B0.',q);
    end
  end
end
