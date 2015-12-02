function [U,iter] = KPM(A,v,F,n,ht,conv,restart,alg,int)
%Skriv en programdefinosjon her
l = size(A,1);
k = length(F);
if max(abs(v)) == 0 || max(max(abs(F))) == 0
    U = sparse(l,k);
    iter = 0;
    return
end
U = zeros(l,k);
iter = 1;
h = norm(v,2);
[Vn,Hn,vnext,hnext] = alg(A,v,n,conv);

[Zn] = int(Hn,[h*F(1,:);sparse(length(Hn)-1,k)],k,ht);

ns = Vn*Zn;
U = U + ns;
diff = hnext;
if restart
    while diff > conv 
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = alg(A,v,n,conv);
        [Zn] = int(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],k,ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
    end
end
end