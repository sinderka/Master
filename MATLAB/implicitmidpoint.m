function [ U ] = implicitmidpoint( A,F,k,ht )
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
mat = inv(speye(n)-A*ht);

U(:,2) = U(:,1) + ht*( A*( U(:,1) + ht/2*( A*U(:,1) + F(:,1) ) ) + 1/2*( F(:,2) + F(:,1) ) );

for i = 3:k
   U(:,i) = ( speye(n) - A*ht )\(U(:,i-2) + 2*ht * ( A*U(:,i-2)/2 + F(:,i-1) )); 
end

end

