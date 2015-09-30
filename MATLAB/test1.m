%function [ output_args ] = untitled( input_args )
%function [ utdata ] = wawesolver( m,n,k,prob,solmeth,conv,para )
% Case: Få 1D sak til å fungere med forskjellige "startverdi
% betingelser"!
clear
close all
m = 100;
k = 100;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k); 
ht = T(2)-T(1);
A = -1/hs^2*gallery('tridiag', m-2);

ant = 1;



% INIT U
U = zeros(m,k);
% INIT F
%F = pi^2*sin(pi*X)'*ones(1,k);
%F = (2-X.*(1-X))'*sin(T);%
%F = -2*ones(m,k);
F = (1+pi^2)*sin(pi*X)'*exp(-T);
% INIT V
%V = zeros(m,1);
%V = (2-X.*(1-X))';
V = -sin(pi*X)';



% Integrate
U(:,1) = sin(pi*X);
%U(:,1) = X.*(X-1);
%U(:,1) = zeros(m,1);
U(2:end-1,2) = U(2:end-1,1) + ht*V(2:end-1) + 0.5*ht^2*(A*U(2:end-1,1) + F(2:end-1,1));

for i = 3:k
    U(2:end-1,i) = 2*U(2:end-1,i-1) - U(2:end-1,i-2) + ht^2*(A*U(2:end-1,i-1) + F(2:end-1,i-1));
end


% Plotting

mesh(U)

%correctsolution = sin(pi*X)'*ones(1,k);
%correctsolution = sin(pi*X)'*cos(pi*T);
%correctsolution = (X.*(X-1))'*ones(1,k);
%correctsolution = (X.*(1-X))'*sin(T);
correctsolution = sin(pi*X)'*exp(-T);

figure(2)
mesh(correctsolution)


figure(3)
mesh(U-correctsolution)
error = max(max(abs(U-correctsolution)))




%end

