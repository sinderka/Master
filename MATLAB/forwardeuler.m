function [ U ] = forwardeuler( A,F,k,ht )
%solves the problem du/dt = Au+F
%indata
% A: mxm matrix
% F: k row
% k: number of points in time
% ht: step size in time
%outdata
% U: the solution
n = size(A,1);
U = zeros(n,k);
for i = 2:k
    U(:,i) = U(:,i-1) + ht*( A*U(:,i-1) + F(:,i-1));
end
end

