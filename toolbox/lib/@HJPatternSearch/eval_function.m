% =================================================================
% Evaluate function
% =================================================================
function F = eval_function(self, x)
  %eval_function This method evaluate the value function and counts the
  %number of evaluations in one iteration

  N   = self.N;
  tol = self.tol;
  lb  = self.lb;
  ub  = self.ub;
            
  % Check bounds and evaluate F to Inf if x is out of bounds
  out_of_lb = any( x < lb+tol );
  out_of_ub = any( ub-tol < x );

  if (out_of_lb || out_of_ub) %&& norm(self.h,inf) < 1e4*self.tol %prod(self.lb <= x) && prod(x <= self.ub)
    F = Inf;
  else
    F = self.fun(x);
  end

  % update statistic
  self.stats.fun_evaluation = self.stats.fun_evaluation + 1;
end
