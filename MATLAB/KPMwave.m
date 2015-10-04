function [U,iter] = KPMwave( Zn,A,V,F,v,k,m,ht,n,conv )
                           %( U,A,V,F,v,k,m,ht,(m-2)^2,conv );

U = zeros(m^2,k);
U(:,1) = Zn(:,1);
Zn = zeros(n,k);
helpvector = zeros((m-2)^2,1);
for qq = 0:m-3
    helpvector(qq*(m-2)+1:qq*(m-2)+m-2) = (qq+1)*m+2:m-1 +(1+ qq)*m;
end

%m = length(A);
%if n == (m-2)^2
    %%%%%Projection method for the full krylov space
%     hn = norm(v,2);
%     [Vm,Hm,~] = Arnoldi(A,v,n,conv);
%     
%     Zn(helpvector,1) = Vm(:,1:n)'*U(helpvector,1);
%     %vector = zeros(n,k); vector(1,:) = Zn(end,:);
%     %[Zn] = integrate(H,vector,n,k,ht,hn);
%     %Zn = doubleintegrate(Zn,V,F,A,ht,n,k );
%     [Zn] = locintegrate( Zn,V,F,Hm,ht,m,k );
%     U(helpvector,:) = Vm(:,1:n)*Zn(helpvector,:);
%     iter = 1;
%else
    %%%%%Projection method for any krylov space
    %U = zeros(m^2,k);
    hn = norm(v,2);
    iter = 1;
    [Vm,Hm,hm] = Arnoldi(A,v,n,conv);
    Zn(:,1) = Vm(:,1:n)'*U(helpvector,1);
    [Zn] = locintegrate( Zn,zeros(n,1),zeros(n,k),Hm,ht,k );
    U(helpvector,2:end) = U(helpvector,2:end) + Vm(:,1:n)*Zn(:,2:end);
    hn = hm; v = Vm(:,end);
    while hn > conv
        [Vm,Hm,hm] = Arnoldi(A,v,n,conv);
        [Zn] = locintegrate( Zn,zeros(n,1),zeros(n,k),Hm,ht,k );
        U(helpvector,2:end) = U(helpvector,2:end) + Vm(:,1:end-1)*Zn(:,2:end);
        hn = hm; v = Vm(:,end);
        iter = iter+1;
        hn
    end
%end

%end
function [Zn] = locintegrate( Zn,V,F,H,ht,k )
    
Zn(:,2) = Zn(:,1) - ht * V(:) + 0.5*ht^2*(H*Zn(:,1) + F(:,1));
for i = 3:k
    Zn(:,i) = 2*Zn(:,i-1) - Zn(:,i-2) + ht^2*(H*Zn(:,i-1) + F(:,i-1));
end

end

end