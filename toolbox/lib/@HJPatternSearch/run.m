% =================================================================
% Run
% =================================================================

function x_sol = run(self)
  %RUN This method run the whole Hooke-Jeeves algorithm. Search
  %is repeated until it fails, then the scal h is reduced. When h
  %is less than a threshold, the method returns the solution.
  
  self.stats.fun_evaluation      = 0; % set to zero the function evaluation counter of each iteration
  self.stats.fun_grad_evaluation = 0; % set to zero the function gradient evaluation counter of each iteration
  
  while norm(self.h,inf) > self.tol && ...
        self.stats.fun_evaluation < self.max_fun_evaluation %k = 1:length(h_vector)
    self.search();
    % If iteration limit is reached,stop. Otherwise, reduce
    % mesh size and continue
    if self.iteration_count >= self.max_iter
      break;
    else
      self.h = self.h*self.rho; % reduce the scale
    end
  end

  if self.verbose
    self.print_info();
  end

  x_sol = self.x;   
end