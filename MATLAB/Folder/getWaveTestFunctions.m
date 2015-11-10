function [Ustart, V,F,correctsolution] = getWaveTestFunctions( prob,X,T )
%Takes 3 agruments;
% prob: a number corresponding to a test problem
% X: a set of points in spacial direction
% T: a set of points in time
%Returns 3 matrices
% F: A list of vectors depending on time
% V: A list of vectors to generate the Krylov space
% correctsolution: The correct solution.
m = length(X); k = length(T);


vec = helpvector(m);

if prob == 1
    sol = @(t,x,y) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);
    u0 = @(x,y) sin(pi*x)*sin(pi*y);
    v0 = @(x,y) 0;
    v1 = @(x,y) 0; f1 = @(t) 0;
    v2 = @(x,y) 0; f2 = @(t) 0;
    V = [getV(v1,X),getV(v2,X)];
    F = [getTime(f1,T);getTime(f2,T)];
elseif prob == 2
    sol = @(t,x,y) (x-1)*x*(y-1)*y*(t^2-t+1);
    u0 = @(x,y) (x-1)*x*(y-1)*y;
    v0 = @(x,y) -(x-1)*x*(y-1)*y;
    v1 = @(x,y) (x-1)*x*(y-1)*y; f1 = @(t) 1;
    v2 = @(x,y) 2*(y*(y-1)+x*(x-1)); f2 = @(t) -(t^2-t+1);
    V = [getV(v1,X),getV(v2,X)];
    F = [getTime(f1,T);getTime(f2,T)];
elseif prob == 3
    sol = @(t,x,y) sin(pi*x)*y*(y-1)*(t^2+1);
    u0 = @(x,y) sin(pi*x)*y*(y-1);
    v0 = @(x,y) 0;
    v1 = @(x,y) 2*sin(pi*x)*y*(y-1); f1 = @(t) 1;
    v2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); f2 = @(t) -(t^2+1);
    V = [getV(v1,X),getV(v2,X)];
    F = [getTime(f1,T);getTime(f2,T)];
elseif prob == 4
    sol = @(t,x,y) sin(pi*x)*y*(y-1)*sin(t);
    u0 = @(x,y) 0;
    v0 = @(x,y) sin(pi*x)*y*(y-1);
    v1 = @(x,y) sin(pi*x)*y*(y-1); f1 = @(t) -sin(t);
    v2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); f2 = @(t) -sin(t);
    V = [getV(v1,X),getV(v2,X)];
    F = [getTime(f1,T);getTime(f2,T)];
elseif prob == 5
    sol = @(t,x,y) sin(pi*x)*sin(pi*y)*sin(sqrt(2)*pi*t);
    u0 = @(x,y) 0;
    v0 = @(x,y) sqrt(2)*pi*sin(pi*x)*sin(pi*y);
    v1 = @(x,y) 0; f1 = @(t) 0;
    v2 = @(x,y) 0; f2 = @(t) 0;
    V = [getV(v1,X),getV(v2,X)];
    F = [getTime(f1,T);getTime(f2,T)];
elseif prob == 6
    sol = @(t,x,y) t^2*x*y*(x-1)*(y-1);
    u0 = @(x,y) 0;
    v0 = @(x,y) 0;
    v1 = @(x,y) x*y*(x-1)*(y-1); f1 = @(t) 2;
    v2 = @(x,y) (x*(x-1)+y*(y-1)); f2 = @(t) -2*t^2;
    V = [getV(v1,X),getV(v2,X)];
    F = [getTime(f1,T);getTime(f2,T)];
elseif prob == 7 % Sjekk!
    sol = @(t,x,y) exp(t*x*y)*(1-x)*x*sin(pi*y);
    u0 = @(x,y) (1-x)*x*sin(pi*y);
    v0 = @(x,y) x^2*y*(1-x)*sin(pi*y);
    f = @(t,x,y)x^2*y^2*exp(x*y*t)*(1-x)*x*sin(pi*y) - ( ...
        t^2*y^2*exp(x*y*t)*(1-x)*x*sin(pi*y) +...
        exp(x*t*y)*sin(pi*y)*(-2) + ...
        2*t*y*sin(pi*y)*(1-2*x) + ...
        t^2*x^2*exp(x*t*y)*(1-x)*x*sin(pi*y) - ...
        pi^2*exp(x*t*y)*sin(pi*y)*(1-x)*x + ...
        2*pi*t*x*exp(x*t*y)*(1-x)*x*cos(pi*y));
    V = speye(2*(m-2)^2);
    F = getSolution(f,X,T);
    F = [sparse((m-2)^2,k);F(vec,:)];
end

Ustart = getInitial(u0,v0,X);
V =  [Ustart,V]; F = [ones(1,k);F];
correctsolution = getSolution(sol,X,T);


function correctsolution = getSolution(sol,X,T)
correctsolution = zeros(m^2,k);
for j = 1:k
    for i = 1:m
        for l = 1:m
            correctsolution(l+(i-1)*m,j) = sol(T(j),X(i),X(l));
        end
    end
end
end
function V =  getInitial(u0,v0,X)
V0 = zeros(m^2,1);
U0 = zeros(m^2,1);
for i = 1:m
    for j = 1:m
        U0(i+(j-1)*m,1) = u0(X(j),X(i));
        V0(i+(j-1)*m,1) = v0(X(j),X(i));
    end
end
V = [U0(vec);V0(vec)];
end
function F = getTime(f,T)
F = zeros(1,k);
for l = 1:k
    F(l) = f(T(l));
end
end
function V =  getV(v,X)
V0 = zeros(m^2,1);
for i = 1:m
    for j = 1:m
        V0(i+(j-1)*m,1) = v(X(j),X(i));
    end
end
V = [sparse((m-2)^2,1);V0(vec)];
end
end