function [U,iter] = KPMwave( Zn,A,V,F,v,k,m,ht,n,conv )

U = zeros(m^2,k);
U(:,1) = Zn(:,1);
Zn = zeros(n,k);
helpvector = zeros((m-2)^2,1);
for qq = 0:m-3
    helpvector(qq*(m-2)+1:qq*(m-2)+m-2) = (qq+1)*m+2:m-1 +(1+ qq)*m;
end

    hn = norm(v,2);
    iter = 1;
    [Vm,Hm,hm] = Arnoldi(A,v,n,conv);
    Zn(:,1) = Vm(:,1:n)'*U(helpvector,1);
    %Ustart =  Vm(:,1:n)*Zn(:,1);
    Vtr = Vm(:,1:n)' *V(helpvector); Ftr = Vm(:,1:n)'*F(helpvector,:);
    [Zn] = locintegrate( Zn,Vtr,Ftr,Hm,ht,k );
    Ustart = Vm(:,1:n)*Zn(:,1);
    U(helpvector,2:end) = Vm(:,1:n)*Zn(:,2:end);
    diff = hm;
    while diff > conv
        hn = hm; v = Vm(:,end);
        [Vm,Hm,hm] = Arnoldi(A,v,n,conv);
        
        vector = spalloc(n,k,k); vector(end,:) = Zn(end,:);
        Zn = zeros(n,k);
        Zn(:,1) = Vm(:,1:n)'*(U(helpvector,1)-Ustart);
        Ustart = Ustart + Vm(:,1:n)*Zn(:,1);
        %Ftr = Vm(:,1:n)'*F(helpvector,:);
        %vector = spalloc(k);
        
        %%% ARG %%% 
        % Ustart blit aldri lik U, hvorfor?
        % Gjøre ting analytiske i morgen!!
        
        
        [Zn] = locintegrate(Zn,sparse(n,1),hm*vector,Hm,ht,k );
        ns =  Vm(:,1:n)*Zn(:,2:end);
        %Ustart = Ustart + Vm(:,1:n)*Zn(:,1);
        diff = max(max(abs(ns)));
        U(helpvector,2:end) = U(helpvector,2:end) + ns;
        iter = iter+1;
        max(U(helpvector,1)-Ustart)
    end



end
function [Zn] = locintegrate( Zn,V,F,H,ht,k )
    
Zn(:,2) = Zn(:,1) - ht * V(:) + 0.5*ht^2*(H*Zn(:,1) + F(:,1));
for i = 3:k
    Zn(:,i) = 2*Zn(:,i-1) - Zn(:,i-2) + ht^2*(H*Zn(:,i-1) + F(:,i-1));
end
end