function [V,H,v,hn]=Arnoldi(A,v,n,conv)
%%%%SKRIV EN PROGRAMDEFINISJON
m = length(A);
n = min(n,m);
V = zeros(m,n+1);
H = zeros(n,n);
V(:,1) = v/(norm(v,2));
for j = 1:1:n
    z = A*V(:,j);
    for i = 1:1:j
        H(i,j) = V(:,i)'*z;
        z = z-H(i,j)*V(:,i);
    end
    hn = norm(z,2);
    if hn < conv
        v = zeros(m,1);
        V = V(:,1:j);
        hn = 0;
        H = H(1:j,1:j);
        return;
    elseif (j+1<=n)
        H(j+1,i) = hn;
        V(:,j+1) = z/(hn);
    else 
        V(:,j+1) = z/(hn);
    end
end
v = V(:,end);
V = V(:,1:n);


end