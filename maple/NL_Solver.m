classdef NL_Solver < handle

  properties (SetAccess = private, Hidden = true)
    tolF2;
    tolGrad;
    iter;
    max_iter;
    ierr;
    f_eval;
    j_eval; 
    more_thuente;
    bfgs;
    trace;
    F0;
    J0;
  end

  methods
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function self = NL_Solver()
      self.tolF2        = 1e-30;
      self.tolGrad      = 1e-12;
      self.max_iter     = 400;
      self.more_thuente = NL_MoreThuente();
      self.bfgs         = NL_BFGS();
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function delete( ~ )
      %% Destroy the C++ class instance
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function set_max_iter( self, miter )
      self.max_iter = miter;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [f,Df] = eval_Phi0( self, x, dir )
      %
      % m(t)  = (1/2) F(x0+t*d)^T W F(x0+t*d) 
      % m'(t) = F(x0+t*d)^T W J(x0+t*d)*d
      %
      % d = -J0^(-1) * F(x0)
      % m'(0) = - F0^T * W J0 * J0^(-1) * F0 = - F0^T * W * F0
      %
      % d = -H * J0^T * W F(x0)
      % m'(0) = - F0^T * W * J0 * H * J0^T * W * F0
      %
      F  = feval( self.f_eval, x );
      Jd = feval( self.j_eval, x )*dir;
      f  = dot(F,F)/2;
      Df = dot(F,Jd);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dir = dir0a( self, J0, F0 )
      n      = length(F0);
      lambda = min(max(max(abs(J0))),max(max(abs(F0))));
      dir = -[J0;lambda*eye(n)]\[F0;zeros(n,1)];
      dir = dir(1:n);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dir = dir0b( self, J0, F0 )
      dir = -self.bfgs.apply(J0.'*F0);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dir = dir0c( self, J0, F0 )
      s = -J0\F0;
      d = -self.bfgs.apply(J0.'*F0);
      if ~all(isfinite(s))
        dir = d;
        return;
      end
      sd = dot(s,d); % = dot(F0,F0)
      s2 = dot(s,s);
      d2 = dot(d,d);
      epsi = 0.1;
      C  = sd^2-epsi^2*d2*s2;
      if C >= 0
        dir = s;
      else
        A = (1-epsi^2)*d2*(d2-2*sd)+C;
        B = 2*((1-epsi^2)*d2*sd-C);
        r = roots([A,B,C]);
        beta = r(1);
        if beta < 0 || beta > 1
          beta = r(2);
        end
        fprintf('BETA %g\n',beta);
        dir = beta*d+(1-beta)*s;
      end
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [f,Df] = eval_Phi1( self, x, dir )
      %
      % m(t)  = (1/2) F(x0+t*d)^T H F(x0+t*d) 
      % m'(t) = F(x0+t*d)^T H J(x0+t*d)*d
      %
      % d     = -J0^(-1) * F0
      % m'(0) = - F0^T H J0*J0^(-1) * F0 = - F0^T H F0
      %
      % d     = -J0^T * H * F0
      % m'(0) = - F0^T H J0*J0^T * H * F0
      %
      F   = feval( self.f_eval, x );
      Jd  = feval( self.j_eval, x ) * dir;
      HF  = self.bfgs.apply(F);
      f   = dot(HF,F)/2;
      Df  = dot(HF,Jd);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dir = dir1a( self, J0, F0 )
      dir = -J0\F0;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dir = dir1b( self, J0, F0 )
      dir = -J0.'*self.bfgs.apply(F0);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [f,Df] = eval_Phi2( self, x, dir )
      %
      % m(t)  = (1/2) F(x0+t*d)^T J0^(-T) * H * J0^(-1) F(x0+t*d) 
      % m'(t) = F(x0+t*d)^T (J0^(-T) * H * J0^(-1)) J(x0+t*d)*d
      %
      % d     = -J0^(-1) * F0
      % m'(0) = - F0^T (J0^(-T) * H * J0^(-1)) J0 * J0^(-1) * F0
      %       = - F0^T (J0^(-T) * H * J0^(-1)) * F0
      %
      F     = feval( self.f_eval, x );
      Jd    = feval( self.j_eval, x ) * dir;
      JJd   = self.J0\Jd;
      JF    = self.J0\F;
      HJF   = self.bfgs.apply(JF);
      f     = dot(HJF,JF)/2;
      Df    = dot(HJF,JJd);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dir = dir2( self, J0, F0 )
      dir = -J0\F0;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [f,Df] = eval_Phi3( self, x, dir )
      %
      % m(t)  = (1/2) F(x0+t*d)^T * H * F(x0+t*d) 
      % m'(t) = F(x0+t*d)^T * H * J(x0+t*d)*d
      %
      % d = -J0^T * H * F0
      % m'(0) = - F0^T * H * J0 * J0^T * H * F0
      %
      F     = feval( self.f_eval, x );
      Jd    = feval( self.j_eval, x ) * dir;
      JJd   = self.J0.'*Jd;
      JF    = self.J0.'*F;
      HJF   = self.bfgs.apply(JF);
      f     = dot(HJF,JF)/2;
      Df    = dot(HJF,JJd);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    function dir = dir3( self, J0, F0 )
      dir = -J0.'*self.bfgs.apply(F0);
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function [x0,ierr] = solve( self, xinit, f_eval, j_eval, ck_eval, do_trace )

      self.iter   = 0;
      self.f_eval = f_eval;
      self.j_eval = j_eval;

      x0 = xinit;
      self.bfgs.setup( length(xinit), 20 );
      alpha = 1;

      self.F0 = feval( self.f_eval, x0 );
      self.J0 = feval( self.j_eval, x0 );
      f0      = dot(self.F0,self.F0)/2;
      g0      = self.J0.'*self.F0;
      self.bfgs.init( x0, g0 );

      if do_trace
        self.trace.iter = 0;
        self.trace.x    = x0;
        self.trace.f2   = f0;
        self.trace.g    = norm(g0,inf);
      end

      for k=1:self.max_iter
        self.iter = k;

        % evaluate merit function and gradient
        if ~feval( ck_eval, x0 )
          fprintf(1, 'NL_CG_Solver: Bad point!\n' );
          ierr = 2;
          return;
        end

        if ~all(isfinite(f0)) %| ~all(isfinite(J0))
          disp(f0);
          fprintf(1, 'NL_CG_Solver: Bad f(x)!\n' );
          ierr = 2;
          return;
        end

        method = 0;
        switch method
        case 0
          d0 = self.dir0a( self.J0, self.F0 ); % 16, NO
          %d0 = self.dir0b( self.J0, self.F0 ); % 16, NO
          %d0 = self.dir0c( self.J0, self.F0 ); % 16, NO
        case 1
          d0 = self.dir1a( self.J0, self.F0 ); % 16, NO
          %d0 = self.dir1b( self.J0, self.F0 ); % 16, NO
        case 2
          d0 = self.dir2( self.J0, self.F0 ); % 16, NO
        case 3
          d0 = self.dir3( self.J0, self.F0 ); % 16, NO
        end

        if ~all(isfinite(d0))
          fprintf(1, 'NL_CG_Solver: bad direction d0!\n' );
          ierr = 10;
          return;
        end

        switch method
        case 0
          fcn        = @(t) self.eval_Phi0( x0+t*d0, d0 );
          [ff0,Dff0] = self.eval_Phi1( x0, d0 );
        case 1
          fcn        = @(t) self.eval_Phi1( x0+t*d0, d0 );
          [ff0,Dff0] = self.eval_Phi1( x0, d0 );
        case 2
          fcn        = @(t) self.eval_Phi2( x0+t*d0, d0 );
          [ff0,Dff0] = self.eval_Phi2( x0, d0 );
        case 3
          fcn        = @(t) self.eval_Phi3( x0+t*d0, d0 );
          [ff0,Dff0] = self.eval_Phi3( x0, d0 );
        end

        if Dff0 >= 0
          fprintf('\nNL_CG_Solver: Df0 = %g!\n',Dff0);
          ierr = 11;
          return;
        end
        stp = 1;
        [step,info,nfev] = self.more_thuente.linesearch(fcn,stp,ff0,Dff0);
        alpha = step.t;

        if ~isfinite(alpha)
          ierr = 3;
          return;
        end

        x0      = x0+alpha*d0;
        self.F0 = feval( self.f_eval, x0 );
        self.J0 = feval( self.j_eval, x0 );
        f0      = dot(self.F0,self.F0)/2;
        g0      = self.J0.'*self.F0;

        normi_f = norm( f0, inf );
        normi_g = norm( g0, inf );

        if do_trace
          self.trace.iter = [self.trace.iter self.iter];
          self.trace.x    = [self.trace.x    x0];
          self.trace.f2   = [self.trace.f2   f0];
          self.trace.g    = [self.trace.g    norm(g0,inf)];
        end

        % check if finished
        if normi_g < self.tolGrad || normi_f < self.tolF2
          fprintf('NL_CG_Solver: iter %3d |f| = %8g |g| = %8g (CONVERGED)\n', ...
                  k, normi_f, normi_g);
          ierr = 0;
          return;
        else
          fprintf('NL_CG_Solver: iter %3d |f| = %8.3g |g| = %8.3g alpha = %8.3g [%4d,%4d]\n',...
                  k, normi_f, normi_g, alpha, info, self.more_thuente.get_iter() );            
        end

        if ~feval( ck_eval, x0 )
          ierr = 7;
          return;
        end

      end
      fprintf(1, 'NL_CG_Solver: too much iterations!\n' );
      ierr = 1;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function it = iterations( self )
      it = self.iter;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
    function trace = get_trace( self )
      trace = self.trace;
    end
    %
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    %
  end
end
