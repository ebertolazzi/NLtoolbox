addpath('minFunc');
addpath('minFunc/compiled');

%x0  = [-0.6;-0.4];
x0  = [-1;0];
x   = linspace(-1,1);
y   = linspace(-1,1);
FUN = Fun1();

[X,Y] = meshgrid(x,y);
XY        = zeros(2,size(X,1),size(X,2));
XY(1,:,:) = X;
XY(2,:,:) = Y;
Z = FUN.fun(XY);

options = {};
options.Display = 'iter';
options.optTol  = 1e-10; % Termination tolerance on the first-order optimality (1e-5)
options.progTol = 1e-10; % Termination tolerance on progress in terms of function/parameter changes (1e-9)

options.Method  = 'qnewton';
% 'sd', 'csd', 'bb', 'cg', 'scg', 'pcg', 'lbfgs',
% 'newton0', 'pnewton0', 'qnewton', 'mnewton'
options.LS_type   = 1;
options.LS_interp = 2;
options.LS_init   = 2;

LST = {'sd', 'csd', 'bb', 'cg', 'pcg', 'lbfgs','qnewton'};
for kkk=1:length(LST)
  subplot(3,3,kkk);
  hold off
  contour(X,Y,log(1+log(1+min(Z,100))),100);
  axis equal
  hold on

  fprintf('\n\n\n%s\n\n\n',LST{kkk});
  options.Method = LST{kkk};
  [x,f,exitflag,output] = minFunc(@funzione,x0,options,FUN);
  plot(output.trace.x(1,:),output.trace.x(2,:),'o-','LineWidth',3);
  title(LST{kkk});
end

function [f,g] = funzione( x, FUN )
  f = FUN.fun(x);
  g = FUN.grad(x);
end