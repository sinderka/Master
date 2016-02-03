function [U,iter,energy1,energy2,energy3] = PM(A,v,F,n,ht,conv,restart,int,figvar,PMint,PMalg)
%Indata
% A: mxm matrix
% v: m vector
% F: k row of timedependant function
% n: real number 0<n<=m
% ht: stepsize in time
% conv: convergence criterion
% restart: A boolean value
% alg: an ortogonalisation algorithm (Arnoldi or SLM)
% int: an integration method (trapezoidal rule)
%outdata
% U: Solution to problem du/dt = Au+v*F
% iter: number of restarts preformed
l = size(A,1);
k = length(F);
if max(abs(v)) == 0 || max(max(abs(F))) == 0
    U = sparse(l,k);
    iter = 0;
    energy1 = 0;
    energy2 = 0;
    energy3 = 0;
    return
end
U = zeros(l,k);
iter = 1;
h = norm(v,2);
[Vn,Hn,vnext,hnext] = PMalg(A,v,n,conv);

if PMint == 1
    Zn = int(Hn,[h*F(1,:);sparse(length(Hn)-1,k)],ht);
elseif PMint == 2
    Zn = expintegrate(Hn,Hn\[h;sparse(length(Hn)-1,1)],0:ht:ht*(k-1));
elseif PMint == 3
    Zn = real(myexpm(full(Hn),Hn\[h;sparse(length(Hn)-1,1)],0:ht:ht*(k-1)));
end

if isequal(PMalg, @Arnoldi)
    energy1 = 0;
    energy2 = 0;
    energy3 = 0;
else
    energy1 = max(energyBIG(A,Zn,v,h,ht,figvar,int));
    energy2 = max(energySMALL(Hn,Vn,Zn,v,h,ht,figvar,int));
    energy3 = max(abs(energyBIG(A,Zn,v,h,ht,0,int)-energySMALL(Hn,Vn,Zn,v,h,ht,0,int)));
end


ns = Vn*Zn;
U = U + ns;
%U1 = int(A,v*F,ht);
diff = hnext;
if restart && diff > conv
    
    
    h = hnext; v = vnext;
    [Vn,Hn,vnext,hnext] = PMalg(A,v,n,conv);
    [Zn] = int(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],ht);
    ns =  Vn*Zn;
    diff = max(max(abs(ns)));
    U = U + ns;
    iter = iter+1;
    
    if isequal(PMalg, @Arnoldi)
        energy1 = 0;
        energy2 = 0;
        energy3 = 0;
    else
        energy1 = max(energyBIG(A,Zn,v,h,ht,figvar,int));
        energy2 = max(energySMALL(Hn,Vn,Zn,v,h,ht,figvar,int));
        energy3 = max(abs(energyBIG(A,Zn,v,h,ht,0,int)-energySMALL(Hn,Vn,Zn,v,h,ht,0,int)));
    end
    
    while diff > conv
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = PMalg(A,v,n,conv);
        [Zn] = int(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
    end
end

end