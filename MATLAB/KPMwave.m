function [U,iter] = KPMwave( Zn,A,V,F,v,k,m,ht,n,conv )
                           %( U,A,V,F,v,k,m,ht,(m-2)^2,conv );

U = zeros(m^2,k);
U(:,1) = Zn(:,1);

helpvector = zeros((m-2)^2,1);
for qq = 0:m-3
    helpvector(qq*(m-2)+1:qq*(m-2)+m-2) = (qq+1)*m+2:m-1 +(1+ qq)*m;
end

%m = length(A);
if n == (m-2)^2
    %%%%%Projection method for the full krylov space
    hn = norm(v,2);
    [Vm,H,~] = Arnoldi(A,v,n);
    
    Zn(helpvector,1) = Vm(:,1:n)'*U(helpvector,1);
    %vector = zeros(n,k); vector(1,:) = Zn(end,:);
    %[Zn] = integrate(H,vector,n,k,ht,hn);
    %Zn = doubleintegrate(Zn,V,F,A,ht,n,k );
    [Zn] = locintegrate( Zn,V,F,H,ht,m,k );
    U(helpvector,:) = Vm(:,1:n)*Zn(helpvector,:);
    iter = 1;
else
    %%%%%Projection method for any krylov space
    U = zeros(length(A),k);
    diff =conv+1;
    hn = norm(v,2);
    iter = 0;
    while diff > conv %hn > conv
        [V,H,hm] = Arnoldi(A,v,n);
        vector = zeros(n,k); vector(1,:) = Zn(end,:);
        %[Zn] = integrate(H,vector,n,k,ht,hn);
        [Zn] = locintegrate( Zn,V,F,A,ht,m,k );
        %(H,F,n,k,ht,hn)
        hn = hm;
        ns = V(:,1:n)*Zn;
        U = U + ns; %V(:,1:n)*Zn;
        v = V(:,end);
        diff = max(max(abs(ns)));
        iter = iter+1;
    end
end

%end
function [Zn] = locintegrate( Zn,V,F,H,ht,m,k )
Zn(helpvector,2) = Zn(helpvector,1) - ht * V(helpvector) + 0.5*ht^2*(H*Zn(helpvector,1) + F(helpvector,1));
for i = 3:k
    Zn(helpvector,i) = 2*Zn(helpvector,i-1) - Zn(helpvector,i-2) + ht^2*(H*Zn(helpvector,i-1) + F(helpvector,i-1));
end

end

end