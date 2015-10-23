function [U,iter] = KPM(A,v,F,k,n,h,ht,conv,restart)
%function [U,iter] = KPMwave(A,U0,G,m,n,k,ht,conv,restarted)

%help = helpvector(m);

if max(F) == min(F)
    vtilde = v*F(1);
else
    vtilde = v*F;
end

U = zeros(h,k);

hn = norm(v,2);
iter = 1;
[Vm,Hm,hm] = Arnoldi(A,v,n,conv);

temp = Vm(:,1:n)'*vtilde;

[Zn] = integrate(Hm,temp,n,k,ht);

ns = Vm(:,1:n)*Zn;
U = U + ns;
diff = hm;
if restart
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

