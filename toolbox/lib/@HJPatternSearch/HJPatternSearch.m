classdef HJPatternSearch < handle
  %% HJPatternSearch class decription
  %
  %  This class implements the Hooke-Jeeves Pattern Search Algorithm.
  %  The references can be found in:
  %
  %  R. Hooke, T. A. Jeeves,
  %  "Direct Search" Solution of Numerical and
  %  Statistical Problems,
  %  Westinghouse Research Laboratories,
  %  Pittsburg, Pennsylvania
  %  
  %  M Bell, Malcolm Pike,
  %  Remark on Algorithm 178: Direct Search,
  %  Communications of the ACM,
  %  Volume 9, Number 9, September 1966, page 684.
  %  
  %  Arthur Kaupe,
  %  Algorithm 178:
  %  Direct Search,
  %  Communications of the ACM,
  %  Volume 6, Number 6, June 1963, page 313.
  %  
  %  FK Tomlin, LB Smith,
  %  Remark on Algorithm 178: Direct Search,
  %  Communications of the ACM,
  %  Volume 12, Number 11, November 1969, page 637-638.
    
  %% Properties - all public
  properties
    rho                % stencile step decreasing factor (must be 0 < rho < 1)
    h                  % scale of the stencil
    fun                % handle to the value function to minimize
    sf                 % stencil failure flag - used to shrink h, sf = 1 means failure
    X0                 % starting point (guess)
    Vmat               % direction matrix for exploration
    x                  % current iteration
    xb                 % base point
    xC                 % stencil center
    fb                 % best value function (fun evaluated in xbc)
    N                  % problem dimension (number of variables)
    d                  % search direction
    tol                % tolerance on the scale h
    lb                 % x lower bound
    ub                 % x upper bound
    max_fun_evaluation % max number of function evaluations
    max_iter           % max number of iterations
    iteration_count    % explore iteration counter
    stats              % fun_evaluation:      function evaluation counter
                       % fun_grad_evaluation: function gradient evaluation counter
    gradient_flag      % if value function gradient is provided this flag is 1
    verbose            % flag to activate info printing
    search_sign        % vector to keep in memory the direction of function value descent from the previous iteration in each direction j 
  end
  
  %% Methods - all public
  methods
    % =================================================================
    % Class constructor
    % =================================================================    
    function self = HJPatternSearch(fun, X0, varargin)
      %HJPatternSearch The constructor initialize the solver
      % parameters and check the inputs when the class is instanciated.
      
      % Switch to overload the constructor
      switch nargin
      case 0
        error('HJPatternSearch needs at lest two inputs arguments');
      case 1
        error('HJPatternSearch needs at lest two inputs arguments');
      case 2
        self.lb = -Inf*ones(size(X0));
        self.ub =  Inf*ones(size(X0));
      case 4
        self.lb = varargin{1};
        self.ub = varargin{2};
      otherwise
        error('Too many input arguments');
      end

      % Initialize parameters and class attributes (properties)
      if isa(fun,'function_handle')
        self.fun = fun; % store function handle
      else
        error('Provided function is not a function handle')
      end
      
      % Check if the provided value function returns also the gradient
      if nargout(self.fun) > 1
        self.gradient_flag = true;
      else
        self.gradient_flag = false;
      end
      
      self.N                              = length(X0);
      self.rho                            = 0.9;           % Initialize stencile step decreasing factor (must be 0 < rho < 1)
      self.X0                             = X0;            % store starting point
      self.Vmat                           = eye(self.N);   % select direction matrix as the identity matrix
      self.x                              = X0;            % initialize current iteration to guess for the first iteration
      self.xC                             = X0;            % initialize stencile center to guess for the first iteration
      self.xb                             = X0;            % initialize base point to guess for the first iteration
      self.sf                             = true;          % initialize stencile failure flag
      self.tol                            = 1e-15;         % initialize scale tolerance
      self.iteration_count                = 0;             % initialize explore iteration counter
      self.stats.fun_evaluation           = 0;             % initialize function evaluation counter
      self.stats.fun_grad_evaluation      = 0;             % initialize function gradient evaluation counter
      self.max_fun_evaluation             = 100000;        % Initialize max number of function evaluations
      self.max_iter                       = 10000;         % Initialize max number of iterations
      self.verbose                        = false;         % Initialize flag to activate info printing
      self.search_sign                    = ones(1,self.N); % Initialize search verse vector to all ones (first verse will be positive for each direction)
      self.fb                             = Inf;
 
      % initialize scale
      for i = 1:self.N
        if self.X0(i) == 0
          self.h(i) = self.rho;
        else
          self.h(i) = self.rho * abs ( self.X0(i) );
        end
      end
    end

    % =================================================================
    % set verbose
    % =================================================================
    function set_verbose(self, tf)
      self.verbose = tf;
    end
    
    % =================================================================
    % info
    % =================================================================
    function print_info(self)
      fprintf('Optimization stats:\n');
      if self.iteration_count >= self.max_iter
        fprintf(...
          'iteration number reached  Max Iteration Limit = %g\n', ...
          self.max_iter ...
        );
      elseif self.stats.fun_evaluation >= self.max_fun_evaluation
        fprintf(...
          'function evaluations %d exceeded the maximum limit [%d]\n', ...
          self.stats.fun_evaluation, self.max_fun_evaluation ...
        );
      else
        fprintf(...
          'mesh size h = %g less than tolerance = %g\n', ...
          norm(self.h,inf), self.tol ...
        );
      end
      fprintf( ...
        [ 'Iterations:                          %d\n' ...
          'Total function evaluations:          %d\n' ...
          'Total function gradient evaluations: %d\n' ], ...
        self.iteration_count, ...
        self.stats.fun_evaluation, ...
        self.stats.fun_grad_evaluation ...
      );
    end

    % =================================================================
    % Evaluate function
    % =================================================================
    F = eval_function(self, x)

    % =================================================================
    % Evaluate function gradient
    % =================================================================
    G = eval_gradient(self, x)

    % =================================================================
    % Explore
    % =================================================================
    explore(self)

    % =================================================================
    % Search
    % =================================================================
    search(self)

    % =================================================================
    % Run
    % =================================================================
    x_sol = run(self)
        
    % =================================================================
    % NRMG Algorithm
    % =================================================================
    R = fnAR(~,X)
  end
end
