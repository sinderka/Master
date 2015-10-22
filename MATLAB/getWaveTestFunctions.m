function [ U0,V0,F1,F2,G1,G2,correctsolution] = getWaveTestFunctions( prob,m,k,X,T )
if prob == 1
    sol = @(t,x,y) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);
    U_0 =@(x,y) sin(pi*x)*sin(pi*y);
    % INIT U
    U0 = zeros(m^2,1);
    for i = 1:m
        for j = 1:m
            U0((i)+(j-1)*m,1) = U_0(X(i),X(j));
        end
    end
    % INIT F
    F1 = sparse(m^2,1);
    F2 = sparse(m^2,1);
    % INIT G
    G1 = sparse(1,k);
    G2 = sparse(1,k);
    % INIT V
    V0 = sparse(m^2,1);
    correctsolution = zeros(m^2,k);
    for j = 1:k
        for i = 1:m
            for l = 1:m
                correctsolution(l+(i-1)*m,j) = sol(T(j),X(i),X(l));
            end
        end
    end
elseif prob == 2
    sol = @(t,x,y) (x-1)*x*(y-1)*y*(t^2-t+1);
    U_0 = @(x,y) (x-1)*x*(y-1)*y;
    V_0 = @(x,y) -(x-1)*x*(y-1)*y;
    F_func1 = @(x,y) 2*x*y*(x-1); F_g1 = @(t) 1;
    F_func2 = @(x,y) (y*(y-1)+x*(x-1)); F_g2 = @(t) 2*(t^2-t+1);
    
    
    V0 = zeros(m^2,1);
    U0 = zeros(m^2,1);
    F1 = zeros(m^2,1);
    F2 = zeros(m^2,1);
    for i = 1:m
        for j = 1:m
            U0((i)+(j-1)*m,1) = U_0(X(i),X(j));
            V0((i)+(j-1)*m,1) = V_0(X(i),X(j));
            F1((i)+(j-1)*m,1) = F_func1(X(i),X(j));
            F2((i)+(j-1)*m,1) = F_func2(X(i),X(j));
        end
    end
    G1 = zeros(1,k);
    G2 = zeros(1,k);
    for l = 1:k
        G1(l) = F_g1(T(l));
        G2(l) = F_g2(T(l));
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
end

