function [ U ] = implicitmidpoint( A,F,k,ht )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

n = size(A,1);
U = zeros(n,k);
mat = inv(speye(n)-A*ht);

%U(:,2) = 2*ht*( F(:,1) ); %U(:,3) = 2*ht*F(:,2); %%1
U(:,2) = U(:,1) + ht*( A*( U(:,1) + ht/2*( A*U(:,1) + F(:,1) ) ) + 1/2*( F(:,2) + F(:,1) ) ); %%2

%U(:,2) = 2*ht*( A*U(:,1) + F(:,1) ); %%3

%U(:,2) = mat*(2*ht*F(:,1));
for i = 3:k
    %U(:,i) = mat*( U(:,i-2) + ht*A*U(:,i-2) + 2*ht*F(:,i-1) ); %%1
   U(:,i) = ( speye(n) - A*ht )\(U(:,i-2) + 2*ht * ( A*U(:,i-2)/2 + F(:,i-1) )); %%2
    %U(:,i) = U(:,i-2) + 2*ht*( A*U(:,i-1) + F(:,i-1) ); %%3
end

end

