function ane()
close all

%%% Test functions
f = @(x) 1/(1+x^2); a = -5; b = 5;corretestimate = 2*atan(5);
%f = @(x) x^2; a = 0; b = 1;corretestimate = 1/3;
%f = @(x) sin(x); a = 0; b = pi;corretestimate = 2;


n = 13; m = 1;

% n = 14;
NewtonCotes(a,b,f,n) % Test problem 1
GaussQuadrature(a,b,f,n) % Test problem 2
alg = @  NewtonCotes; compositeIntegration(a,b,f,n,m,alg) % Test problem 4
alg = @  GaussQuadrature; compositeIntegration(a,b,f,n,m,alg) % Test problem 4

% Problem 3
plotresult(a,b,f,corretestimate)

% Problem 4
m = 5; compositPlotresult(a,b,f,m,corretestimate);
end
function su = compositeIntegration(a,b,f,n,m,alg) % composit function for problem 4.

su = 0;
for i = 0:m-1
    atemp = a + (b-a)/m*i;
    btemp = a + (b-a)/m*(i+1);
    su = su + alg(atemp,btemp,f,n);
end

end

function compositPlotresult(a,b,f,m,corretestimate) % Plot function for problem 4. Note the @ symbols, what do you think they do?
maxn = 13;
N = zeros(1,maxn);
G = zeros(1,maxn);
for i = 1:maxn
    N(i) = abs(compositeIntegration(a,b,f,i,m,@NewtonCotes)-corretestimate);
    G(i) = abs(compositeIntegration(a,b,f,i,m,@GaussQuadrature)-corretestimate);
end
figure(2)
semilogy(1:maxn,N)
hold on
loglog(1:maxn,G)
legend('Composit Newton-Cotes','Composit Gauss-Quadrature')
xlabel('n');
ylabel('error');
end

function plotresult(a,b,f,corretestimate) % Plot results problem 3
maxN = 14;
N = zeros(1,maxN);
G = zeros(1,maxN);
for i = 1:maxN
    N(i) = abs(NewtonCotes(a,b,f,i)-corretestimate);
    G(i) = abs(GaussQuadrature(a,b,f,i)-corretestimate);
end
figure(1)
semilogy(1:maxN,N)
hold on
loglog(1:maxN,G)
legend('Newton-Cotes','Gauss-Quadrature')
xlabel('n');
ylabel('error');
end

function su = NewtonCotes(a,b,f,n) % Problem 1
% https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas
% Does not always converge
load('NewtonCotes.mat');
weights = (b-a)*weights(n,:);
h = (b-a)/n;
su = 0;
for i = 0:n
    su = su + weights(i+1)*f(a+h*i);
end
end
function su = GaussQuadrature(a,b,f,n) % Problem 2
% https://en.wikipedia.org/wiki/Gaussian_quadrature
% Seams to always work fine!
load('GaussQuadrature.mat')
weights = (b-a)/2*weights(n,:);
points = (b-a)/2*points(n,:)+(a+b)/2;
su = 0;
for i = 0:n
    su = su + weights(i+1)*f(points(i+1));
end
end