%function [ utdata ] = wawesolver( m,n,k,prob,solmeth,conv,para )
clear
close all
m = 120;
k = 120;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k); 
ht = T(2)-T(1);
A = -1/hs^2*gallery('tridiag', m-2);

ant = 1;


U = zeros(m,k);
U(:,1) = sin(pi*X);
U(2:end-1,2) = U(2:end-1,1) + 0.5*ht^2*(A*U(2:end-1,1));

for i = 3:k
    U(2:end-1,i) = 2*U(2:end-1,i-1) - U(2:end-1,i-2) + ht^2*(A*U(2:end-1,i-1));

end

mesh(U)

correctsolution = sin(pi*X)'*cos(pi*T);


figure(2)
mesh(correctsolution)


figure(3)
mesh(U-correctsolution)
error = max(max(abs(U-correctsolution)))


