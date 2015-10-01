function [ U ] = doubleintegrate( U,V,F,A,ht,m,k )

helpvector = zeros((m-2)^2,1);
for i = 0:m-3
    helpvector(i*(m-2)+1:i*(m-2)+m-2) = (i+1)*m+2:m-1 +(1+ i)*m;
end


U(helpvector,2) = U(helpvector,1) - ht * V(helpvector) + 0.5*ht^2*(A*U(helpvector,1) + F(helpvector,1));
for i = 3:k
    U(helpvector,i) = 2*U(helpvector,i-1) - U(helpvector,i-2) + ht^2*(A*U(helpvector,i-1) + F(helpvector,i-1));
end

end

