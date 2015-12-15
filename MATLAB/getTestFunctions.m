function [Ustart, V,F,correctsolution] = getTestFunctions( prob,X,T,eqn )
%Takes 4 agruments;
% prob: a number corresponding to a test problem
% X: a set of points in spacial direction
% T: a set of points in time
% eqn: spesifies which equation we are solving
%Returns 4 matrices
% Ustart: The initial value of the testproblem
% F: A list of vectors depending on time
% V: A list of vectors to generate the Krylov space
% correctsolution: The correct solution.

%%% TODO
% Legg til maxwell problemer
m = length(X); k = length(T);


[vec,~] = helpvector(m,eqn);
if strcmp(eqn,'heat') % Sjekk!!!
    if prob == 1
        sol = @(t,x,y)  t/(t+1)*x*(x-1)*y*(y-1);
        u0 = @(x,y) 0;
        f1 = @(t) 1/(t+1)^2;   v1 = @(x,y) x*(x-1)*y*(y-1);
        f2 = @(t) -t/(t+1);     v2 = @(x,y) 2*x*(x-1)+2*y*(y-1);
        V = [getV(v1,X),getV(v2,X)];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 2
        sol = @(t,x,y) exp(x*y)*y*(y-1)*sin(pi*x)*t*cos(t);
        u0 = @(x,y) 0;
        f1 = @(t) cos(t)-t*sin(t);            v1 = @(x,y) (y-1)*y*exp(x*y)*sin(pi*x);
        f2 = @(t) -t*cos(t);           v2 = @(x,y) exp(x*y)*(x^2*(y-1)*y+x*(4*y-2)+2)*sin(pi*x) +(y-1)*y^3*exp(x*y)*sin(pi*x)+...
            2*pi*(y-1)*y^2*exp(x*y)*cos(pi*x)-pi^2*(y-1)*y*exp(x*y)*sin(pi*x);
        V = [getV(v1,X),getV(v2,X)];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 3
        sol = @(t,x,y) (t^2+1)*(x-1)*x*y*(y-1);
        u0 = @(x,y) (x-1)*x*y*(y-1);
        f1 = @(t) 2*t;   v1 = @(x,y) (x-1)*x*y*(y-1);
        f2 = @(t) -(t^2+1);     v2 = @(x,y) 2*(x-1)*x+2*(y-1)*y;
        V = [getV(v1,X),getV(v2,X)];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 4
        sol = @(t,x,y) (x-1)*x*sin(pi*y)*(t^2+1);
        u0 = @(x,y) (x-1)*x*sin(pi*y);
        f1 = @(t) 2*t;   v1 = @(x,y) (x-1)*x*sin(pi*y);
        f2 = @(t) -(t^2+1);     v2 = @(x,y) 2*sin(pi*y) - pi^2*(x-1)*x*sin(pi*y);
        V = [getV(v1,X),getV(v2,X)];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 5
        sol = @(t,x,y)sin(x*y*t)*(x-1)*(y-1);
        f = @(t,x,y) t^2*(x-1)*(y-1)*y^2*sin(t*x*y)-2*t*(y-1)*y*cos(t*x*y)+(x-1)*x*(y-1)*y*cos(t*x*y)-t*(x-1)*x*(2*cos(t*x*y)-t*x*(y-1)*sin(t*x*y));
        V = speye((m-2)^2);
        F = getSolution(f,X,T);
        F = F(vec,:);
    end
    Ustart = getInitial(u0,X);
    V =  [Ustart,V]; F = [ones(1,k);F];
    correctsolution = getSolution(sol,X,T);
elseif strcmp(eqn,'wave')
    if prob == 1
        sol = @(t,x,y) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);
        u0 = @(x,y) sin(pi*x)*sin(pi*y);
        v0 = @(x,y) 0;
        v1 = @(x,y) 0; f1 = @(t) 0;
        v2 = @(x,y) 0; f2 = @(t) 0;
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 2
        sol = @(t,x,y) (x-1)*x*(y-1)*y*(t^2-t+1);
        u0 = @(x,y) (x-1)*x*(y-1)*y;
        v0 = @(x,y) -(x-1)*x*(y-1)*y;
        v1 = @(x,y) (x-1)*x*(y-1)*y; f1 = @(t) 1;
        v2 = @(x,y) 2*(y*(y-1)+x*(x-1)); f2 = @(t) -(t^2-t+1);
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 3
        sol = @(t,x,y) sin(pi*x)*y*(y-1)*(t^2+1);
        u0 = @(x,y) sin(pi*x)*y*(y-1);
        v0 = @(x,y) 0;
        v1 = @(x,y) 2*sin(pi*x)*y*(y-1); f1 = @(t) 1;
        v2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); f2 = @(t) -(t^2+1);
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 4
        sol = @(t,x,y) sin(pi*x)*y*(y-1)*sin(t);
        u0 = @(x,y) 0;
        v0 = @(x,y) sin(pi*x)*y*(y-1);
        v1 = @(x,y) sin(pi*x)*y*(y-1); f1 = @(t) -sin(t);
        v2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); f2 = @(t) -sin(t);
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 5
        sol = @(t,x,y) sin(pi*x)*sin(pi*y)*sin(sqrt(2)*pi*t);
        u0 = @(x,y) 0;
        v0 = @(x,y) sqrt(2)*pi*sin(pi*x)*sin(pi*y);
        v1 = @(x,y) 0; f1 = @(t) 0;
        v2 = @(x,y) 0; f2 = @(t) 0;
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 6
        sol = @(t,x,y) t^2*x*y*(x-1)*(y-1);
        u0 = @(x,y) 0;
        v0 = @(x,y) 0;
        v1 = @(x,y) x*y*(x-1)*(y-1); f1 = @(t) 2;
        v2 = @(x,y) (x*(x-1)+y*(y-1)); f2 = @(t) -2*t^2;
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 7
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
    
    Ustart = [getInitial(u0,X);getInitial(v0,X)];
    V =  [Ustart,V]; F = [ones(1,k);F];
    correctsolution = getSolution(sol,X,T);
elseif strcmp(eqn,'maxwell1D')
    if prob == 1
        u0 = @(x) 0;
        v0 = @(x) cos(pi*x);
        
        %v0 = @(x) sin(pi*x);
        %u0 = @(x) cos(pi*x);
        solE = @(t,x) sin(pi * x) * sin(pi * t);
        %solB = @(t,x) cos(pi * x) * cos(pi * t);
        F = ones(1,k);
        %V = zeros(2*(m-2),1);
        %V(m-2) = 1/(2*hs); V(1) = -1/(2*hs);
        %F(2,:) = cos(pi*T);
        %V(m-1) = -1/(2*hs); V(end) = -1/(2*hs);
        correctsolution = zeros(m,k);
        for j = 1:k
            for i = 1:m
                correctsolution(i,j) = solE(T(j),X(i));
                %correctsolution(m+i,j) = solB(T(j),X(i));
            end
        end
        %     elseif prob == 2
        %         u0 = @(x) exp(-100*(x-0.5)^2);
        %         v0 = @(x) exp(-100*(x-0.5)^2);
        %         correctsolution = sparse(m,k);
        %         F = ones(1,k);
        %     elseif prob == 3
        %         u0 = @(x) cos(x);
        %         v0 = @(x) cos(x);
        %         solE = @(t,x) cos(x-t);
        %         F = ones(1,k);
        %         correctsolution = zeros(m,k);
        %         for j = 1:k
        %             for i = 1:m
        %                 correctsolution(i,j) = solE(T(j),X(i));
        %                 %correctsolution(m+i,j) = solB(T(j),X(i));
        %             end
        %         end
    end
    U0 = zeros(m-2,1); V0 = zeros(m,1);
    %V0(1) = V0(T(1)); V0(end) = v0(T(end));
    for i = 1:m
        V0(i) = v0(X(i));
    end
    
    for i = vec
        U0(i-1) = u0(X(i));
    end
    Ustart = [U0;V0];
    V = [Ustart];
    
elseif strcmp(eqn,'maxwell3D')
    %fyll inn!
    
elseif strcmp(eqn,'random')
    Ustart = rand(2*(m-2)^2,1);
    V = Ustart;
    F = ones(1,k);
    correctsolution = sparse(m^2,k);
elseif strcmp(eqn,'semirandom')
    if prob == 1
        try
            load('semirandomV.mat','V');
        catch
            V = -1;
        end
        
        if size(V,1) ~= 2*(m-2)^2
            V = rand(2*(m-2)^2,1);
            save('semirandomV.mat','V');
        end
        Ustart = V;
        F = ones(1,k);
        correctsolution = sparse(m^2,k);
    elseif prob == 2
        try
            load('semirandomV.mat','V');
        catch
            V = -1;
        end
        
        if size(V,1) ~= 2*(m-2)^2
            V = rand(2*(m-2)^2,1);
            save('semirandomV.mat','V');
        end
        
        try 
            load('semirandomF.mat','F')
        catch
            F = -1;
        end
        Ustart = V;
        if size(F,2) ~= k
            F = rand(1,k);
            save('semirandomF.mat','F');
        end

        correctsolution = sparse(m^2,k);
    end
    
end

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
    function U0 =  getInitial(u0,X)
        U0 = zeros(m^2,1);
        for i = 1:m
            for j = 1:m
                U0(i+(j-1)*m,1) = u0(X(j),X(i));
            end
        end
        U0 = U0(vec);
    end
    function F = getTime(f,T)
        F = zeros(1,k);
        for l = 1:k
            F(l) = f(T(l));
        end
    end
    function V0 =  getV(v,X)
        V0 = zeros(m^2,1);
        for i = 1:m
            for j = 1:m
                V0(i+(j-1)*m,1) = v(X(j),X(i));
            end
        end
        V0 = V0(vec);
    end
end