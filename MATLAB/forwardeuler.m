function [ U ] = forwardeuler( A,F,k,ht )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = size(A,1);
U = zeros(n,k);
for i = 2:k
    U(:,i) = U(:,i-1) + ht*( A*U(:,i-1) + F(:,i-1));
end
end

