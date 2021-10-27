% =================================================================
% Evaluate function gradient
% =================================================================
function G = eval_gradient(self, x)
  %eval_gradient This method evaluate the gradient of the value function and counts the
  %number of evaluations in one iteration

  % Check bounds and evaluate G to Inf if x is out of bounds
  if all(self.lb <= x) && all(x <= self.ub)
    [~, G] = self.fun(x);
  else
    G = Inf*ones(self.N,1);
  end

  % update statistic
  self.stats.fun_grad_evaluation = self.stats.fun_grad_evaluation + 1;
end
