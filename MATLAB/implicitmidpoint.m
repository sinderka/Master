function [ U ] = implicitmidpoint( A,F,ht )
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
mat = inv(speye(n)-A*ht);

U(:,2) = mat*( U(:,1) + ht*F(:,2) );
for i = 3:2:k-1
    U(:,i) = U(:,i-1) + ht*( A*U(:,i-1) + F(:,i-1) );
    U(:,i+1) = mat*(U(:,i) + ht*F(:,i+1));
end

end

