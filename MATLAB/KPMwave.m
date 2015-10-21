function [U,iter] = KPMwave2(A,U0,m,n,k,ht,conv)

v = helpvector(m);
%v2 = v1+m^2;

U0tilde = A*U0;
%Residual = U0tr;
U = zeros(2*(m-2)^2,k);
%U(helpvector,:) = Zn(helpvector,1)*ones(1,k);
%Zn = zeros(n,k);

hn = norm(v,2);
iter = 1;
%temp = ones(162,1);% temp(1) = 1;
[Vm,Hm,hm] = Arnoldi(A,U0tilde,n,conv);
[Vm1,Hm1] = krylov(A,U0tilde,n);
%U0tilde = (J*U0);
%temp = spalloc(n,k,k); 
%temp = Vm(:,1:n)'*(U0tilde*ones(1,k)); 
%temp = [temp(1,:);sparse(2*(m-2)^2-1,k)];

%temp = [Vm(:,1)'*U0tilde*ones(1,k);sparse(n-1,k)];
temp = Vm(:,1:n)'*U0tilde*ones(1,k);
%temp2 = Vm(:,1:n)'*J*Vm(:,1:n);

[Zn] = locintegrate2(Hm,temp,n,k,ht);

ns = Vm(:,1:n)*Zn;
U = U + ns;
diff = hm;
while diff > conv
    hn = hm; v = Vm(:,end);
    [Vm,Hm,hm] = Arnoldi(A,v,n,conv);

    temp = [hn*Zn;sparse(n-1,k)];
    
    %temp2 = Vm(:,1:n)'*J*Vm(:,1:n);
    [Zn] = locintegrate2(Hm,temp,n,k,ht);
    ns =  Vm(:,1:n)*Zn;
    
    diff = max(max(abs(ns)));
    U = U + ns;
    iter = iter+1;
    
end



end


function Zn = locintegrate2(H,F,n,k,ht)
%%%%%Solves the equation z'-H*z=F Numerically
%%% ht is stepsize, n is the size of the square matrix H
%%% k is the number of time-steps.

    mat = inv(eye(n)-ht/2*H);
    Zn = zeros(n,k);
    Zn(:,2) = mat*(ht/2*(F(:,1)+F(:,2)));
    for j = 3:1:k
        Zn(:,j) = mat*(Zn(:,j-1)+ht/2*(H*Zn(:,j-1)+(F(:,j)+F(:,j-1))));
    end
end