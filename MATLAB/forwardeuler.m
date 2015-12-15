function [ U ] = forwardeuler( A,F,ht )
%solves the problem du/dt = A*u+F
%indata
% A: mxm matrix
% F: mxk matrix
% ht: step size in time
%outdata
% U: the solution
k = size(F,2);
n = size(A,1);
U = zeros(n,k);
for i = 2:k
    U(:,i) = U(:,i-1) + ht*( A*U(:,i-1) + F(:,i-1));
end
end

