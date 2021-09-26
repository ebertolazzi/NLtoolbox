classdef NL_StrongWolfe < handle

  properties (SetAccess = private, Hidden = true)
    c_1_vec;
    c_2_vec;
    rho;
    tot_iter;
    max_iterations;
    phi_0;
    Dphi_0;
    A_list;
    B_list;
  end

  methods
 
    function self = NL_StrongWolfe()
      self.c_1_vec        = [0.1,0.01,1e-4]; % Armijo condition
      self.c_2_vec        = [0.1,0.5,0.9];   % second (strong) Wolfe condition
      self.rho            = 2;               % bracket growth
      self.max_iterations = [5,10,20];
    end
  
    % a_lo = a_{i - 1}
    % a_hi = a_{i}
    function res = interpolate( self, a_lo, a_hi, phi_lo, phi_hi, Dphi_lo, Dphi_hi )
      d1  = Dphi_lo + Dphi_hi - 3 * (phi_lo - phi_hi) / (a_lo - a_hi);
      d2  = sqrt(d1 * d1 - Dphi_lo * Dphi_hi);
      tmp = (Dphi_hi + d2 - d1) / (Dphi_hi - Dphi_lo + 2 * d2);
      tmp = max(0.1,min(0.9,tmp));
      res = a_hi - (a_hi - a_lo) * tmp;
    end

    function res = zoom( self, a_lo_in, a_hi_in, phi, Dphi )
      a_lo = a_lo_in; phi_lo = feval( phi, a_lo); Dphi_lo = feval( Dphi, a_lo );
      a_hi = a_hi_in; phi_hi = feval( phi, a_hi); Dphi_hi = feval( Dphi, a_hi );

      % Shrink bracket
      for kkk=1:3
        c1_Dphi_0 = self.c_1_vec(kkk)*self.Dphi_0;
        c2_Dphi_0 = self.c_2_vec(kkk)*self.Dphi_0;
        for iteration=1:self.max_iterations(kkk)
          self.tot_iter = self.tot_iter + 1;

          self.A_list = [ self.A_list, a_lo ];
          self.B_list = [ self.B_list, a_hi ];

          % Interpolate a_j
          if a_lo < a_hi
            a_j = self.interpolate( a_lo, a_hi, phi_lo, phi_hi, Dphi_lo, Dphi_hi );
          else
            % TODO: Check if this is needed
            a_j = self.interpolate( a_hi, a_lo, phi_hi, phi_lo, Dphi_hi, Dphi_lo );
          end

          % Evaluate phi(a_j)
          phi_j  = feval( phi,  a_j );
          Dphi_j = feval( Dphi, a_j );

          % Check Armijo
          if (phi_j > self.phi_0 + a_j * c1_Dphi_0) || (phi_j > phi_lo)
            a_hi    = a_j;
            phi_hi  = phi_j;
            Dphi_hi = Dphi_j;
          else

            if abs(Dphi_j) <= -c2_Dphi_0
              res = a_j;
              return;
            end

            if Dphi_j * (a_hi - a_lo) >= 0
              a_hi    = a_lo;
              phi_hi  = phi_lo;
              Dphi_hi = Dphi_lo; 
            end

            a_lo    = a_j;
            phi_lo  = phi_j;
            Dphi_lo = Dphi_j; 
          end
        end
      end

      % Quasi-error response
      res = a_j;
    end

    % `StrongWolfe`: This linesearch algorithm guarantees that 
    % the step length satisfies the (strong) Wolfe conditions.
    % See Nocedal and Wright - Algorithms 3.5 and 3.6
    % This algorithm is mostly of theoretical interest,
    % users should most likely
    % `MoreThuente`, `HagerZhang` or `BackTracking`.

    function res = linesearch( self, phi, Dphi, alpha0 )

      self.tot_iter = 0;

      self.phi_0  = feval( phi, 0 );
      self.Dphi_0 = feval( Dphi, 0 );

      self.A_list = [];
      self.B_list = [];

      % Step-sizes
      a_0     = 0;
      a_im1   = a_0;
      a_i     = alpha0;
      phi_im1 = self.phi_0;

      for kkk=1:3
        c1_Dphi_0 = self.c_1_vec(kkk)*self.Dphi_0;
        c2_Dphi_0 = self.c_2_vec(kkk)*self.Dphi_0;
        for iteration=1:self.max_iterations(kkk)

          self.tot_iter = self.tot_iter + 1;

          phi_i = feval( phi, a_i );

          % Test Wolfe conditions
          if (phi_i > self.phi_0 + a_i * c1_Dphi_0) || ( (phi_i >= phi_im1) && ( a_im1 > 0 ) )
%             ( (phi_i >= phi_im1) && ( (iteration > 1) || (kkk > 1) ) )
            a_star = self.zoom( a_im1, a_i, phi, Dphi );
            res = [a_star, feval( phi, a_star )];
            return;
          end

          Dphi_i = feval( Dphi, a_i );

          % Check condition 2
          if abs(Dphi_i) <= -c2_Dphi_0
            res = [a_i, phi_i];
            return;
          end

          % Check condition 3
          if Dphi_i >= 0 % FIXME untested!
            a_star = self.zoom( a_i, a_im1, phi, Dphi );
            res = [a_star, feval( phi, a_star )];
            return;
          end

          self.A_list = [ self.A_list, a_im1 ];
          self.B_list = [ self.B_list, a_i ];

          % Choose a_iplus1 from the interval (a_i, a_max)
          a_im1 = a_i;
          a_i   = a_i * self.rho;

          % Update phi_im1
          phi_im1 = phi_i;
        end
      end

      % Quasi-error response TODO make this error instead
      printf("\n\nlinesearch failed!\n\n");
      res = [a_im1, feval( phi, a_im1 )];
    end

    function res = get_iter( self )
      res = self.tot_iter;
    end

    function doplot( self, f_in, df_in, alpha, amax, do_interval )
      f0  = feval( f_in, 0 );
      df0 = feval( df_in, 0 );
      x   = 0:amax/400:amax;
      y   = f_in(x);
      plot( x, y );
      hold on;
      plot([0,amax],[f0,f0+amax*self.rho*df0],'Color','LimeGreen');
      plot([0,amax],[f0,f0+amax*self.sigma*df0],'Color','blue');
      plot(alpha,f_in(alpha),'or');
      if do_interval
        for kkk=1:length(A_list)
          plot([A_list(kkk),B_list(kkk)],[-kkk/5,-kkk/5],'Color','black');
        end
      end
    end
  end
end
