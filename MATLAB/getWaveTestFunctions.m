function [ U,V,F,correctsolution] = getWaveTestFunctions( prob,m,k,X,T )
if prob == 1
    sol = @(t,x,y) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);
    U_0 =@(x,y) sin(pi*x)*sin(pi*y);
    % INIT U
    U = zeros(m^2,k);
    for i = 1:m
        for j = 1:m
            U((i)+(j-1)*m,1) = U_0(X(i),X(j));
        end
    end
    % INIT F
    F = zeros(m^2,k);
    % INIT V
    V = zeros(m^2,1);
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

