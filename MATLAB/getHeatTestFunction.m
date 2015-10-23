function [U0,F1,F2,G1,G2,correctsolution] = getHeatTestFunction(prob,m,k,X,T)
func1 = 0; func2 = 0; func3 = 0; func4 = 0; sol = 0;

if prob == 1
    sol = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1);
    u0 = @(x,y) 0;
    func1 = @(t) 1/(t+1).^2;   func2 = @(x,y) x*(x-1)*y.*(y-1);
    func3 = @(t) -t/(t+1);     func4 = @(x,y) 2*x.*(x-1)+2*y.*(y-1);
    
elseif prob == 2
    sol = @(t,x,y) exp(x.*y).*y.*(y-1).*sin(pi*x).*t.*cos(t);
    u0 = @(x,y) 0;
    func1 = @(t) cos(t)-t*sin(t);            func2 = @(x,y) (y-1).*y.*exp(x.*y).*sin(pi*x);
    func3 = @(t) -t.*cos(t);           func4 = @(x,y) exp(x.*y).*(x.^2.*(y-1).*y+x.*(4*y-2)+2).*sin(pi*x) +(y-1).*y.^3.*exp(x.*y).*sin(pi*x)+...
        2*pi*(y-1).*y.^2.*exp(x.*y).*cos(pi*x)-pi^2*(y-1).*y.*exp(x.*y).*sin(pi*x);
    %elseif prob == 3
    %    sol = @(t,x,y)sin(x.*y.*t).*(x-1).*(y-1);
    %    func1 = @(t,x,y) t.^2.*(x-1).*(y-1).*y.^2.*sin(t.*x.*y)-2*t.*(y-1).*y.*cos(t.*x.*y)+(x-1).*x.*(y-1).*y.*cos(t.*x.*y)-t.*(x-1).*x.*(2*cos(t.*x.*y)-t.*x.*(y-1).*sin(t.*x.*y));
elseif prob == 3
    sol = @(t,x,y) (t^2+1)*(x-1)*x*y*(y-1);
    u0 = @(x,y) (x-1)*x*y*(y-1);
    func1 = @(t) 2*t;   func2 = @(x,y) (x-1)*x*y*(y-1);
    func3 = @(t) -(t^2+1);     func4 = @(x,y) 2*(x-1)*x+2*(y-1)*y;
elseif prob == 4
    sol = @(t,x,y) (x-1)*x*sin(pi*y)*(t^2+1);
    u0 = @(x,y) (x-1)*x*sin(pi*y);
    func1 = @(t) 2*t;   func2 = @(x,y) (x-1)*x*sin(pi*y);
    func3 = @(t) -(t^2+1);     func4 = @(x,y) 2*sin(pi*y) - pi^2*(x-1)*x*sin(pi*y);
end

F1 = zeros(m^2,1);
F2 = zeros(m^2,1);
U0 = zeros(m^2,1);
for i = 1:m
    for j = 1:m
        U0(i+(j-1)*m,1) = u0(X(j),X(i));
        F1(i+(j-1)*m,1) = func2(X(j),X(i));
        F2(i+(j-1)*m,1) = func4(X(j),X(i));
    end
end
G1 = zeros(1,k);
G2 = zeros(1,k);
for i = 1:k
    G1(i) = func1(T(i));
    G2(i) = func3(T(i));
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