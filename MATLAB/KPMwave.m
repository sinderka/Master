function [U,iter] = KPMwave(A,U0,G,m,n,k,ht,conv,restarted)

v = helpvector(m);

if max(G) == min(G)
    U0tilde = U0*G(1);
else
    U0tilde = U0*G;
end

U = zeros(2*(m-2)^2,k);

hn = norm(v,2);
iter = 1;
[Vm,Hm,hm] = Arnoldi(A,U0,n,conv);

temp = Vm(:,1:n)'*U0tilde;

[Zn] = integrate(Hm,temp,n,k,ht);

ns = Vm(:,1:n)*Zn;
U = U + ns;
diff = hm;
if restarted
    while diff > conv
        hn = hm; v = Vm(:,end);
        [Vm,Hm,hm] = Arnoldi(A,v,n,conv);
        
        temp = [hn*Zn(end,:);sparse(n-1,k)];
        
        [Zn] = integrate(Hm,temp,n,k,ht);
        ns =  Vm(:,1:n)*Zn;
        
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
        
    end
end


end


% function Zn = locintegrate2(H,F,n,k,ht)
% %%%%%Solves the equation z'-H*z=F Numerically
% %%% ht is stepsize, n is the size of the square matrix H
% %%% k is the number of time-steps.
%
% mat = inv(eye(n)-ht/2*H);
% Zn = zeros(n,k);
% Zn(:,2) = mat*(ht/2*(F(:,1)+F(:,2)));
% for j = 3:1:k
%     Zn(:,j) = mat*(Zn(:,j-1)+ht/2*(H*Zn(:,j-1)+(F(:,j)+F(:,j-1))));
% end
% end