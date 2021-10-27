% =================================================================
% Search
% =================================================================

function search(self)
  %SEARCH This method call the explore method on the first
  %iteration and then continue to call explore until a stencil
  %fails. In the case of a stencile failure, it tries once to go
  %back of half a step along the search direction by setting xC equal to the base point xb.
  %If the stencil fails again, it exits the while loop and sf is set to zero in order
  %to signal that a reduction of h is necessary.
  
  self.iteration_count = self.iteration_count + 1; % augment counter
  % Print info
  if self.verbose
    %line = '----------------------------------------------------------------------------------------------------';
    line = '-------------------------------------------------------------------------';
    fprintf('%s\n',line');
    fprintf( ...
      'Iteration=%-5d f(x)=%-10.5g scale=%-10.5g #f eval=%-5d #G eval=%-5d\n',...
      self.iteration_count, self.fb, norm(self.h,inf), ...
      self.stats.fun_evaluation, self.stats.fun_grad_evaluation ...
    );
    for ii=1:self.N
      fprintf( '%d\t%g\n', ii, self.x(ii) );
    end
    fprintf('%s\n',line');
  end
  
  % Explore the first iteration point
  self.xb = self.x; % Move the base point to the current iteration
  self.xC = self.x; % Set the new stencil center
  self.sf = true;
  
  self.explore();

  keep = true;
  % Continue exploring until stencil failure or exceed of
  while self.sf && keep && self.stats.fun_evaluation < self.max_fun_evaluation
    self.d  = self.x - self.xb; % Compute search direction
    self.xb = self.x;          % Move the base point to the current iteration
    self.xC = self.x + self.d;  % Set the new stencil center at a distance 2d from xb

    self.explore(); % Explore the stencil centered in xC
      
    % If there is a stencil failure, move xC back to the past
    % iteration x (which is equal to xb ) and explore
    if ~self.sf
      self.xC = self.xb;
      self.explore(); % Explore the stencil centered in xC
    end
    
    % If there is a stencil failure again, the search method
    % will terminate (and h must be reduced to continue the algorithm)

    % Check whether the current best point is changed or it is
    % still the initial base point xb %% check total number of iterations
    
    % if the new point is different in at least one of the direction,
    % break and use the flag keep to assert the new point is different
    % from the previous
    if any( 0.5 * self.h < abs(self.x-self.xb) )
      keep = false;
    end
  end
end
