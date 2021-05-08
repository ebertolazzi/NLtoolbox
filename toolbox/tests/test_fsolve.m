clc;
clear all;
close all;

global nl;

nl = NLtest();

%res = nl.listall()
%options = optimoptions(@fsolve,'Display','none','SpecifyObjectiveGradient',true,'PlotFcn',@optimplotfirstorderopt);
algo = {'trust-region-dogleg', 'trust-region','levenberg-marquardt'};


stat = {};
kkk  = 0 ;

for ntest=1:nl.numTest()
  nl.select(ntest) ;
  options = optimoptions( @fsolve,...
    'Display',                  'none',...
    'Algorithm',                algo{3}, ...
    'SpecifyObjectiveGradient', true, ...
    'FunctionTolerance',        1e-10, ...
    'MaxIterations',            100,...
    'OptimalityTolerance',      1e-10, ...
    'CheckGradients',           false,...
    'JacobPattern',             nl.pattern(), ...
    'FiniteDifferenceType',     'central', ...
    'FiniteDifferenceStepSize', eps^(1/3.5) ...
  );
  for ng=1:nl.numGuess()
    
    fprintf('test N.%4d GUESS N.%2d PROBLEM = %s\n',ntest,ng,nl.name());

    x0 = nl.guess(ng);
    [x,F,exitflag,output,JAC] = fsolve(@fun,x0,options);
    kkk = kkk+1 ;
    stat{kkk}.exitflag   = exitflag;
    stat{kkk}.title      = nl.name();
    stat{kkk}.iterations = output.iterations;
    stat{kkk}.x          = x;
    stat{kkk}.nguess     = ng;
    stat{kkk}.n          = nl.numActive();
  end
end

ok  = 0;
qok = 0;
for k=1:kkk
  conv = 'DIVERGE';
  if stat{k}.exitflag == 1
    conv = 'YES    ';
    ok = ok+1;
  end
  if stat{k}.exitflag == 4
    conv = 'NO     ';
    qok = qok+1;
  end

  fprintf('%4d Conv:%7s Test N:%4d Guess:%2d IT:%4d [%s]\n', ...
           k, ...
           conv, ...
           stat{k}.n, ...
           stat{k}.nguess, ...
           stat{k}.iterations, ...
           stat{k}.title);
  % disp(stat{k}.x.');
end

fprintf('ntest = %d\n',kkk);
fprintf('ok   = %d\n',ok);
fprintf('qok  = %d\n',qok);
fprintf('fail = %d\n',kkk-qok-ok);

fprintf('\n\nAll done folks\n\n');

function [F,J] = fun(x)
  global nl ;
  F = nl.evalF(x);
  J = nl.evalJF(x);
end