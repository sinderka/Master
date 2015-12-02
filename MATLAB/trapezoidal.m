function [ U ] = trapezoidal( A,F,k,ht )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n = size(A,1);
U = zeros(n,k);
mat = inv(speye(n)-A*ht/2);
for i = 2:k
    U(:,i) = mat*(U(:,i-1) + ht/2*A*U(:,i-1)+ht/2*(F(:,i)+F(:,i-1)));
end
end

