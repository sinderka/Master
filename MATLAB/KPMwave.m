function [U,iter] = KPMwave( Zn,A,v,k,ht,n,conv )

m = length(A);
if n == m
    %%%%%Projection method for the full krylov space
    hn = norm(v,2);
    [Vm,H,~] = Arnoldi(A,v,n);
    vector = zeros(n,k); vector(1,:) = Zn(end,:);
    %[Zn] = integrate(H,vector,n,k,ht,hn);
    Zn = doubleintegrate(Zn,V,F,A,ht,n,k );
    U = Vm(:,1:n)*Zn;
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
        [Zn] = integrate(H,vector,n,k,ht,hn);
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

