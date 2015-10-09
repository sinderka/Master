function [U,iter] = KPMwave( Zn,A,V,F,v,k,m,ht,n,conv )

helpvector = zeros((m-2)^2,1);
for qq = 0:m-3
    helpvector(qq*(m-2)+1:qq*(m-2)+m-2) = (qq+1)*m+2:m-1 +(1+ qq)*m;
end

U0tr = A*Zn(helpvector,1);
Residual = U0tr;
U = zeros(m^2,k);
U(helpvector,:) = Zn(helpvector,1)*ones(1,k);
Zn = zeros(n,k);

    hn = norm(v,2);
    iter = 1;
    [Vm,Hm,hm] = Arnoldi(A,v,n,conv);
    Vtr = Vm(:,1:n)' *V(helpvector);
    temp = Vm(:,1:n)'*U0tr*ones(1,k);
    [Zn] = locintegrate( Zn,Vtr,temp,Hm,ht,k );
    
    %%% Fungerer nesten, hva mangler? %%%
    
    U(helpvector,:) = U(helpvector,:) + Vm(:,1:n)*Zn;
    diff = hm;
    while diff > conv
        hn = hm; v = Vm(:,end);
        [Vm,Hm,hm] = Arnoldi(A,v,n,conv);
        vector = zeros(n,k); vector(end,:) = Zn(end,:);
        vector = vector + Vm(:,1:n)'*U0tr*ones(1,k);
        Zn = zeros(n,k);
        
        [Zn] = locintegrate(Zn,sparse(n,1),hn*vector,Hm,ht,k );
        ns =  Vm(:,1:n)*Zn;

        diff = max(max(abs(ns)));
        U(helpvector,:) = U(helpvector,:) + ns;
        iter = iter+1;
        
    end



end
function [Zn] = locintegrate( Zn,V,F,H,ht,k )
    
Zn(:,2) = Zn(:,1) - ht * V(:) + 0.5*ht^2*(H*Zn(:,1) + F(:,1));
for i = 3:k
    Zn(:,i) = 2*Zn(:,i-1) - Zn(:,i-2) + ht^2*(H*Zn(:,i-1) + F(:,i-1));
end
end