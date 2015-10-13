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
%temp = spalloc(n,k,k); 
%temp = Vm(:,1:n)'*(U0tilde*ones(1,k)); 
%temp = [temp(1,:);sparse(2*(m-2)^2-1,k)];
temp = [Vm(:,1)'*U0tilde*ones(1,k);sparse(2*(m-2)^2-1,k)];
[Zn] = locintegrate2(Hm,temp,n,k,ht);

%%% Fungerer nesten, hva mangler? %%%
ns = Vm(:,1:n)*Zn;
U = U +ns;
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

% Må nok testes litt
function Zn = locintegrate2(H,F,n,k,ht) %%% Intigrere direkte
%%%%%Solves the equation z'-H*z=hn*e_1*F Numerically
%%% ht is stepsize, n is the size of the square matrix H
%%% k is the number of time-steps.
%%% Zn(:,1) = 0.
    mat = inv(eye(n)-ht/2*H);
    %e1 = zeros(n,1); e1(1) = 1;
    Zn = zeros(n,k);
    Zn(:,2) = mat*(ht/2*(F(:,1)+F(:,2)));
    for j = 3:1:k
        Zn(:,j) = mat*(Zn(:,j-1)+ht/2*(H*Zn(:,j-1)+(F(:,j)+F(:,j-1))));
    end
end