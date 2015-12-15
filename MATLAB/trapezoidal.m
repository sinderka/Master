function [ U ] = trapezoidal( A,F,ht )
%solves the problem du/dt = Au+F
%indata
% A: mxm matrix
% F: mxk matrix
% ht: step size in time
%outdata
% U: the solution
k = size(F,2);
n = size(A,1);
U = zeros(n,k);
mat = inv(speye(n)-A*ht/2);
for i = 2:k
    U(:,i) = mat*(U(:,i-1) + ht/2*A*U(:,i-1)+ht/2*(F(:,i)+F(:,i-1)));
end
end

