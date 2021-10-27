clear all;
close all;
clc;

%% Load test function
global nl;
% Instance test library
nl = NLtest();

%% Cycle throught the tests
num_test                       = 0;
Matlab_patternsearch_conv      = 0;
HJSolverWithGradient_converged = 0;
HJSolver_converged             = 0;

%% Statistics
conv_tolerance = 1e-6; % converged treshold

for ntest = 1:nl.numberOfTests
  nl.select(ntest);
    
  % Select only tests under selected number of variables
  if nl.neq() > 10
    continue;
  end

  % Cycle on all guesses of a single problem
  for ng=1:nl.numGuess()

    fprintf('test N.%4d GUESS N.%2d PROBLEM = %s\n', ntest, ng, nl.name() ); % print test name and guess

    %% Initial guess
    X0       = nl.guess(ng);
    [lb, ub] = nl.bbox();

    %% My solver no gradient
    HJSolver = HJPatternSearch( @myfun, X0, lb, ub );

    %% My solver with gradient
    HJSolverWithGradient = HJPatternSearch( @myfunwithgradient, X0, lb, ub );

    tic
    x_sol = HJSolver.run;
    mysolver_time = toc;

    HJSolver.print_info()

    value_function = myfun(x_sol);
    if value_function < conv_tolerance
      HJSolver_converged = HJSolver_converged+1;
      fprintf('\nHJSolver Converged \n');
    else
      fprintf('\nHJSolver Not Converged \n');
    end

    tic
    x_sol_gradient = HJSolverWithGradient.run();
    mysolver_time_with_gradient = toc;

    HJSolverWithGradient.print_info()

    value_function_with_gradient = myfunwithgradient(x_sol_gradient);
    if value_function_with_gradient < conv_tolerance
      HJSolverWithGradient_converged = HJSolverWithGradient_converged+1;
      fprintf('\nHJSolverWithGradient Converged \n');
    else
      fprintf('\nHJSolverWithGradient Not Converged \n');
    end

    %% Matlab
    tic
    x_sol_matlab = patternsearch(@myfun,X0,[],[],[],[],lb,ub);
    matlab_solver_time = toc;

    value_function_matlab = myfun(x_sol_matlab);        
    if value_function_matlab < conv_tolerance
      Matlab_patternsearch_conv = Matlab_patternsearch_conv+1;
      fprintf('\npatternsearch Converged\n');
    else
      fprintf('\npatternsearch Not Converged\n');
    end
    num_test = num_test+1;
  end
  % pause()
end

%% Print statistics

fprintf('\nNumber of test: %d\n',num_test);

line = '---------------------------------------------------------------------------------';
fprintf( ...
  '%s\nConverged:     HJSolver = %-4d\tHJSolverWithGradient=%-4d\tpatternsearch=%-4d\n', ...
  line, HJSolver_converged, HJSolverWithGradient_converged, Matlab_patternsearch_conv ...
);

fprintf( ...
  '%s\nNot Converged: HJSolver = %-4d\tHJSolverWithGradient=%-4d\tpatternsearch=%-4d\n', ...
  line, ...
  num_test - HJSolver_converged, ...
  num_test - HJSolverWithGradient_converged, ...
  num_test - Matlab_patternsearch_conv ...
);
fprintf( '%s\n', line );

function [NF] = myfun(x)
  global nl ;
  NF = norm(nl.evalF(x));
end

function [NF,G] = myfunwithgradient(x)
  global nl ;
  F  = nl.evalF(x);
  NF = norm(F);
  % Compute the gradient of the 2-norm of multivariate function F(x) as
  % (JF(x)^T*F(x))/norm2(F(x))
  G = (nl.evalJF(x).') * F;
end
