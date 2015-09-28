function [func1,func2,func3,func4,sol] = getTestFunction(eqn,prob)
func1 = 0; func2 = 0; func3 = 0; func4 = 0; sol = 0;

if strcmp(eqn,'heat')
    if prob == 1
        sol = @(t,x,y)  t/(t+1)*x*(x-1)*y.*(y-1);
        func1 = @(t) 1/(t+1).^2;   func2 = @(x,y) x*(x-1)*y.*(y-1);
        func3 = @(t) -t/(t+1);     func4 = @(x,y) 2*x.*(x-1)+2*y.*(y-1);
    elseif prob == 2
        sol = @(t,x,y) exp(x.*y).*y.*(y-1).*sin(pi*x).*t.*cos(t);
        func1 = @(t) cos(t)-t*sin(t);            func2 = @(x,y) (y-1).*y.*exp(x.*y).*sin(pi*x);
        func3 = @(t) -t.*cos(t);           func4 = @(x,y) exp(x.*y).*(x.^2.*(y-1).*y+x.*(4*y-2)+2).*sin(pi*x) +(y-1).*y.^3.*exp(x.*y).*sin(pi*x)+...
            2*pi*(y-1).*y.^2.*exp(x.*y).*cos(pi*x)-pi^2*(y-1).*y.*exp(x.*y).*sin(pi*x);
    elseif prob == 3
        sol = @(t,x,y)sin(x.*y.*t).*(x-1).*(y-1);
        func1 = @(t,x,y) t.^2.*(x-1).*(y-1).*y.^2.*sin(t.*x.*y)-2*t.*(y-1).*y.*cos(t.*x.*y)+(x-1).*x.*(y-1).*y.*cos(t.*x.*y)-t.*(x-1).*x.*(2*cos(t.*x.*y)-t.*x.*(y-1).*sin(t.*x.*y));
    end
    %display('Error')
    
elseif strcmp(eqn,'wave')
    if prob == 1
        func1 = @(t,x) sin(pi*t)*sin(pi*x);
    end
    %display('Error')
elseif strcmp(eqn,'maxwell')
    if prob == 1
        func1 = 0;
    end
    %display('Error')
end
%display('Error')


end