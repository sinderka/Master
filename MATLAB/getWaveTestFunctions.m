function [ U0,V0,F,correctsolution] = getWaveTestFunctions( prob,m,k,X,T )
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
    F = sparse(m^2,k);
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
% elseif prob == 2
%     sol = @(t,x,y) (x-1)*x*(y-1)*y*(t^2-t+1);
%     U_0 = @(x,y) (x-1)*x*(y-1)*y;
%     V_0 = @(x,y) -(x-1)*x*(y-1)*y;
%     F_func = @(t,x,y) 2*x*y*(x-1)+2*(t^2-t+1)*(y*(y-1)+x*(x-1));
%     
%     U0 = 
    
    
end


end

