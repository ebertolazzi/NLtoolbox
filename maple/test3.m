addpath('fun');

%x0  = [1;1;1];
x0  = [0;1e-20;0];
FUN = Potra3();

CG = NL_Solver();

F  = @(x) FUN.F(x(1),x(2),x(3));
JF = @(x) FUN.J(x(1),x(2),x(3));
CK = @(x) true;

[x,ierr] = CG.solve( x0, F, JF, CK, true );

trace = CG.get_trace();
x

figure('Position',[1,1,1200,600]);

subplot(1,2,1)
axis tight;
semilogy(trace.iter,trace.f2,'o-','LineWidth',3);
title('F2');

subplot(1,2,2)
axis tight;
semilogy(trace.iter,trace.g,'o-','LineWidth',3);
title('GRAD');
