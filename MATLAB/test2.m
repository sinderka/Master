%function [ utdata ] = wavesolver( m,n,k,prob,solmeth,conv,para )
% Få til en naiv implementasjon av 2D wave equation
clear
close all
m = 50;
k = 400;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k);
ht = T(2)-T(1);
A = -1/hs^2*gallery('poisson', m-2);

ant = 1;


% INIT U
U = zeros(m^2,k);
for i = 1:m
    for j = 1:m
        U((i)+(j-1)*m,1) = sin(pi*X(i))*sin(pi*X(j));
    end
end
% INIT F
F = zeros(m^2,k);
% INIT V
V = zeros(m^2,1);

helpvector = zeros((m-2)^2,1);
for i = 0:m-3
    helpvector(i*(m-2)+1:i*(m-2)+m-2) = (i+1)*m+2:m-1 +(1+ i)*m;
end


    U(helpvector,2) = U(helpvector,1) - ht * V(helpvector) + 0.5*ht^2*(A*U(helpvector,1) + F(helpvector,1));
    for i = 3:k
        U(helpvector,i) = 2*U(helpvector,i-1) - U(helpvector,i-2) + ht^2*(A*U(helpvector,i-1) + F(helpvector,i-1));
    end
sol = @(t,x,y) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);


error = max(max(max(abs(U-correctsolution))))
video(U,m,k,0.05)
video(correctsolution,m,k,0.05)
video(U-correctsolution,m,k,0.05)

