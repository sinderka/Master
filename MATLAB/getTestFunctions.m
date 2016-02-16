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
% Første vektoren i V skal være initsialbetingelser! Uansett!!!!!!

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
        F = getSolution(X,T,f);
        F = F(vec,:);
    end
    Ustart = getInitial(u0,X);
    V =  [Ustart,V]; F = [ones(1,k);F];
    correctsolution = getSolution(X,T,sol);
elseif strcmp(eqn,'wave') % Kun problem 1 og 2 fungerer
    
    
    if prob == 1
        a1 = 1; b1 = 2; c1 = sqrt(a1^2+b1^2);
        sol1 = @(t,x,y)  sin(a1*pi*x)*sin(b1*pi*y)*cos(c1*pi*t);
        sol2 = @(t,x,y) -sin(a1*pi*x)*sin(b1*pi*y)*c1*pi*sin(c1*pi*t);
        u0 = @(x,y) sin(a1*pi*x)*sin(b1*pi*y);
        
        v0 = @(x,y) 0;
        v1 = @(x,y) 0; f1 = @(t) 0;
        v2 = @(x,y) 0; f2 = @(t) 0;
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 2
        return
        sol1 = @(t,x,y) (x-1)*x*(y-1)*y*(t^2-t+1);
        sol2 = @(t,x,y) (x-1)*x*(y-1)*y*(2*t-1);
        u0 = @(x,y) (x-1)*x*(y-1)*y;
        v0 = @(x,y) -(x-1)*x*(y-1)*y;
        v1 = @(x,y) 2*(x-1)*x*(y-1)*y; f1 = @(t) 1;
        v2 = @(x,y) 2*(y*(y-1)+x*(x-1)); f2 = @(t) -(t^2-t+1);
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 3
        sol1 = @(t,x,y) sin(pi*x)*y*(y-1)*(t^2+1);
        sol2 = @(t,x,y) sin(pi*x)*y*(y-1)*(2*t);
        u0 = @(x,y) sin(pi*x)*y*(y-1);
        v0 = @(x,y) 0;
        v1 = @(x,y) 2*sin(pi*x)*y*(y-1); f1 = @(t) 1;
        v2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); f2 = @(t) -(t^2+1);
        V = [[sparse((m-2)^2,1);getV(v1,X)],[sparse((m-2)^2,1);getV(v2,X)]];
        F = [getTime(f1,T);getTime(f2,T)];
    elseif prob == 4
        sol1 = @(t,x,y) (t-x) + 5*(t+x);
        sol2 = @(t,x,y) 6;
        u0 = 4*x;
        v0 = 0;
        v1 = 0;     f1 = @(t) 0;
        v2 = 0;     f2 = @(t) 0;
        %         sol = @(t,x,y) sin(pi*x)*y*(y-1)*sin(t);
        %         u0 = @(x,y) 0;
        %         v0 = @(x,y) sin(pi*x)*y*(y-1);
        %         v1 = @(x,y) sin(pi*x)*y*(y-1); f1 = @(t) -sin(t);
        %         v2 = @(x,y) sin(pi*x)*(2-pi^2*y*(y-1)); f2 = @(t) -sin(t);
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
    correctsolution = getSolution(X,T,sol1,sol2);
elseif strcmp(eqn,'maxwell1D') % Må sjekkes
    if prob == 1
        u0 = @(x) 0;
        v0 = @(x) cos(pi*x);
        
        %v0 = @(x) sin(pi*x);
        %u0 = @(x) cos(pi*x);
        solE = @(t,x) sin(pi * x) * sin(pi * t);
        solB = @(t,x) cos(pi * x) * cos(pi * t);
        F = [ones(1,k);zeros(1,k)];
        %V = zeros(2*(m-2),1);
        %V(m-2) = 1/(2*hs); V(1) = -1/(2*hs);
        %F(2,:) = cos(pi*T);
        %V(m-1) = -1/(2*hs); V(end) = -1/(2*hs);
        correctsolution = zeros(m,k);
        for j = 1:k
            for i = 1:m
                correctsolution(i,j) = solE(T(j),X(i));
                correctsolution(m+i,j) = solB(T(j),X(i));
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
    
    for i = 2:m-1
        U0(i-1) = u0(X(i));
    end
    Ustart = [U0;V0];
    V = [Ustart,zeros(2*m-2,1)];
    
elseif strcmp(eqn,'maxwell3D') % uferdig
    %fyll inn!
    
elseif strcmp(eqn,'random') % ubrukelig
    Ustart = rand(2*(m-2)^2,1);
    V = Ustart;
    F = ones(1,k);
    correctsolution = sparse(2*m^2,k);
elseif strcmp(eqn,'semirandom') % Kjenner ikke løsningen på problemet den løser
    if prob == 1
        try
            load('semirandomV.mat','V');
        catch
            V = -1;
        end
        
        if size(V,1) ~= 2*(m-2)^2
            V = [rand(2*(m-2)^2,1),zeros(2*(m-2)^2,1)];
            %V = [rand(2*(m-2)^2,1)];
            save('semirandomV.mat','V');
        end
        
        Ustart = V(:,1);
        
        %F = ones(2,k);
        F = ones(2,k);
        correctsolution = sparse(2*m^2,k);
        correctsolution(1,:) = ones(1,k);
    elseif prob == 2
        try
            load('semirandomV2.mat','V');
        catch
            V = -1;
        end
        
        if size(V,1) ~= 2*(m-2)^2
            V = [zeros(2*(m-2)^2,1),rand(2*(m-2)^2,1)];
            save('semirandomV2.mat','V');
        end
        
        try
            load('semirandomF2.mat','F')
        catch
            F = -1;
        end
        Ustart = V(:,1);
        if size(F,2) ~= k
            F = [zeros(1,k);rand(1,k)];
            save('semirandomF.mat','F');
        end
        
        correctsolution = sparse(2*m^2,k);
        correctsolution(1,:) = ones(1,k);
    end
elseif strcmp(eqn,'eigen')
    try
        load('semirandomV.mat','V');
    catch
        V = -1;
    end
    
    if size(V,1) ~= 2*(m-2)^2
        V = rand(2*(m-2)^2,1);
        save('semirandomV.mat','V');
    end
    %Ustart = V;
    F = ones(1,k);
    Ustart = V(:,1);
    correctsolution = zeros(2*m^2,k);
    
    load('eigenQ.mat','Q');
    eigenvals = eigenvalues(2*(m-2)^2-1);
    D = gallery('tridiag',-eigenvals,zeros(2*(m-2)^2,1),eigenvals);
    %[a,b] = eigs(D);
    expmat = exponfunction(D);
    expmat1 = expm(D);
    
    correctsolution1 = zeros(2*m^2,k);
    for i = 1:k
        %correctsolution(vec,i) = Q*expm(D*T(i))*Q'*V(:,1);
        correctsolution(vec,i) =  Q*expmat *Q'*V(:,1);
        correctsolution1(vec,i) = Q*expmat1*Q'*V(:,1);
    end
end

    function correctsolution = getSolution(X,T,sol1,sol2)
        
        if nargin == 4
            correctsolution = zeros(2*m^2,k);
            for j = 1:k
                for i = 1:m
                    for l = 1:m
                        correctsolution(l+(i-1)*m,j) = sol1(T(j),X(i),X(l));
                        correctsolution(m^2+l+(i-1)*m,j)=sol2(T(j),X(i),X(l));
                    end
                end
            end
        else
            correctsolution = zeros(m^2,k);
            for j = 1:k
                for i = 1:m
                    for l = 1:m
                        correctsolution(l+(i-1)*m,j) = sol1(T(j),X(i),X(l));
                    end
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
        if strcmp(eqn,'wave')
            U0 = U0(vec(1:length(vec)/2));
        else
            U0 = U0(vec);
        end
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
        if strcmp(eqn,'wave')
            V0 = V0(vec(1:length(vec)/2)); % denne lager problemer senere!
        else
            V0 = V0(vec);
        end
        
    end
end