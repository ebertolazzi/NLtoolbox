addpath('fun');

x0  = [-0.6;-0.4];
FUN = Fun2();

CG = NL_Solver();

F  = @(x) FUN.F(x(1),x(2));
JF = @(x) FUN.J(x(1),x(2));
CK = @(x) true;

[x,ierr] = CG.solve( x0, F, JF, CK, true );

trace = CG.get_trace();
x

figure('Position',[1,1,800,800]);

subplot(2,2,[1,2])
x   = linspace(-1,1);
y   = linspace(-1,1);
[X,Y] = meshgrid(x,y);
XY    = zeros(2,size(X,1),size(X,2));
Z     = FUN.fun(X,Y);

hold off
ZZ = min(Z-min(min(Z)),100);
contour(X,Y,log(1+ZZ),100);
axis equal
hold on
plot(trace.x(1,:),trace.x(2,:),'o-','LineWidth',3);
axis tight
%title(LST{kkk});

subplot(2,2,3)
semilogy(trace.iter,trace.f2,'o-','LineWidth',3);
axis tight;
title('F2');

subplot(2,2,4)
semilogy(trace.iter,trace.g,'o-','LineWidth',3);
axis tight;
title('GRAD');
