function U = expintegrate(A,U0,T)
%solves the problem du/dt = Au, u(0) = 0
%indata
% A: mxm matrix
% u0: initial condition
% k: number of points in time
% ht: step size in time
%outdata
% U: the solution
% Option 1 (nothing)
%expA = expm(A); %Option 2
%[V,D]= eig(A);%Option 3

U = zeros(length(A),length(T));
for i = 1:length(T)
U(:,i) = expm(A*T(i))*U0-U0;
end