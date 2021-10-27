% =================================================================
% Explore
% =================================================================

function explore(self)
  %EXPLORE This method explore all points on the stencil center at xt = xC and
  % updates the current iteration x to the current best point xcb. If
  % the current best point xcb is worse than the base point xb,
  % the current iteration x will remain constant (x = xb) and
  % stencil failure flag sf will be set to zero.
            
  % Initialize
  self.fb = self.eval_function(self.xb);
  self.d  = 0;
  self.sf = false;
  xcb     = self.xb;
  xt      = self.xC; % temporary point representing the center of the stencil

  % ----------------------------------------------------------------------------------------
  % If the gradient of the value function is provided,
  % align first search direction with gradient
  if self.gradient_flag
    % Evaluate gradient in xt, change sign (decrease direction) and normalize it
    grad = -self.eval_gradient(xt);
    grad = grad./norm(grad);
  
    % Compute n-dimensional rotation matrix from grad to
    % x1 direction
    R = self.fnAR(grad');
      
    % Set direction matrix so that first direction is in gradient direction
    self.Vmat = R';      
  end
  % ----------------------------------------------------------------------------------------
  % Cycle on all stencil directions
    
  for j = 1:size(self.Vmat,2)
    dirh   = self.h(j) * self.Vmat(:,j);
    p      = xt + self.search_sign(j)*dirh;
    fb_tmp = self.eval_function(p);
    if fb_tmp >= self.fb
      p = xt - self.search_sign(j)*dirh; % try the opposite direction
      fb_tmp = self.eval_function(p);
      if ~self.gradient_flag
        self.search_sign(j) = - self.search_sign(j); % change priority of search direction to the opposite verse
      end
    end
    
    % Update temporary and current best point before checking
    % the remaining directions j
    if fb_tmp < self.fb
      xt     = p; % move temporary point
      xcb    = p; % new current best point
      self.fb = fb_tmp; % new best value function
      self.sf = true; % update stencil failure flag  
    end
  end
  if self.sf     
    self.x = xcb; %update the current iteration to current best
    % Increase mesh in case of new best point to speed up
    % search
    % if self.gradient_flag
    %   self.h = self.h*1.01; % Increase the stencile
    % end
  end  
end
