%function [ utdata ] = wawesolver( m,n,k,prob,solmeth,conv,para )
clear
close all
m = 30;
k = 30;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k); 
ht = T(2)-T(1);
%A = -ht^2/hs^2*gallery('tridiag', m-2);

A = -1/hs^2*gallery('tridiag', m-2);
%A = -1/hs^2*gallery('tridiag', m);

%U = 0;
ant = 1;


U = zeros(m,k+1);
U(:,2) = sin(pi*X);
U(:,1) = U(:,2);
%F(1:length(X)) = cos(pi*X);
%invMat = inv(ht^2*A);
c = 0.5;
for i = 3:k+1
    U(2:end-1,i) = 2*U(2:end-1,i-1) - U(2:end-1,i-2) +c*ht^2*(A*U(2:end-1,i-1));
    c = 1;

    %U(:,i) = 2*U(:,i-1) - U(:,i-2) +65*ht^2*(A\U(:,i-1));
    % DET ER EN FEIL I LINJEN OVER!!!!!!!!!!
    %U(2,i) = U(3,i);
    %U(end-1,i) = U(end-2,i);
end
%U = reshape(U,m,k);
U(:,1) = [];
mesh(U)

%f = @(t,x) cos(pi*t)*sin(pi*x);

%correctsolution = zeros(m,k);
%for i = 1:k
%    for j = 2:m-1
%        correctsolution(j,i) = f(T(i),X(j));
%    end
%end
correctsolution = sin(pi*X)'*cos(pi*T);


figure(2)
mesh(correctsolution)


figure(3)
mesh(U/max(max(U))-correctsolution/max(max(correctsolution)))
error = max(max(abs(U-correctsolution)))
minimum = min(min(U))
min = min(U(:,end))
%min = min(U(:,end-1))
%end

