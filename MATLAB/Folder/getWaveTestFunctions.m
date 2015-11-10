function [ U0,V0,F1,F2,G1,G2,correctsolution] = getWaveTestFunctions( prob,X,T )
%Skriv en programdefinosjon her
m = length(X); k = length(T);

if prob == 1
    sol = @(t,x,y) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);
    u0 = @(x,y) sin(pi*x)*sin(pi*y);
    v0 = @(x,y) 0;
    f1 = @(x,y) 0; g1 = @(t) 0;
    f2 = @(x,y) 0; g2 = @(t) 0;
    
elseif prob == 2
    sol = @(t,x,y) (x-1)*x*(y-1)*y*(t^2-t+1);
    u0 = @(x,y) (x-1)*x*(y-1)*y;
    v0 = @(x,y) -(x-1)*x*(y-1)*y;
    f1 = @(x,y) (x-1)*x*(y-1)*y; g1 = @(t) 1;
    f2 = @(x,y) 2*(y*(y-1)+x*(x-1)); g2 = @(t) -(t^2-t+1);
    
elseif prob == 3
    sol = @(t,x,y) sin(pi*x)*y*(y-1)*(t^2+1);
    u0 = @(x,y) sin(pi*x)*y*(y-1);
    v0 = @(x,y) 0;
    f1 = @(x,y) 2*sin(pi*x)*y*(y-1); g1 = @(t) 1;
    f2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); g2 = @(t) -(t^2+1);
elseif prob == 4
    sol = @(t,x,y) sin(pi*x)*y*(y-1)*sin(t);
    u0 = @(x,y) 0;
    v0 = @(x,y) sin(pi*x)*y*(y-1);
    f1 = @(x,y) sin(pi*x)*y*(y-1); g1 = @(t) -sin(t);
    f2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); g2 = @(t) -sin(t);
elseif prob == 5
    sol = @(t,x,y) sin(pi*x)*sin(pi*y)*sin(sqrt(2)*pi*t);
    u0 = @(x,y) 0;
    v0 = @(x,y) sqrt(2)*pi*sin(pi*x)*sin(pi*y);
    f1 = @(x,y) 0; g1 = @(t) 0;
    f2 = @(x,y) 0; g2 = @(t) 0;
elseif prob == 6
    sol = @(t,x,y) t^2*x*y*(x-1)*(y-1);
    u0 = @(x,y) 0;
    v0 = @(x,y) 0;
    f1 = @(x,y) x*y*(x-1)*(y-1); g1 = @(t) 2;
    f2 = @(x,y) (x*(x-1)+y*(y-1)); g2 = @(t) -2*t^2;
end
V0 = zeros(m^2,1);
U0 = zeros(m^2,1);
F1 = zeros(m^2,1);
F2 = zeros(m^2,1);
for i = 1:m
    for j = 1:m
        U0(i+(j-1)*m,1) = u0(X(j),X(i));
        V0(i+(j-1)*m,1) = v0(X(j),X(i));
        F1(i+(j-1)*m,1) = f1(X(j),X(i));
        F2(i+(j-1)*m,1) = f2(X(j),X(i));
    end
end
G1 = zeros(1,k);
G2 = zeros(1,k);
for l = 1:k
    G1(l) = g1(T(l));
    G2(l) = g2(T(l));
end


correctsolution = zeros(m^2,k);
for j = 1:k
    for i = 1:m
        for l = 1:m
            correctsolution(l+(i-1)*m,j) = sol(T(j),X(i),X(l));
        end
    end
end
end
%function
