function [U,iter,energy1,energy2] = SLM(A,v,F,n,ht,conv,restart,int,figvar,PMint)
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
if max(abs(v)) == 0 || max(abs(F)) == 0
    U = sparse(l,k);
    iter = 0;
    energy1 = 0;
    energy2 = 0;
    return
end
U = zeros(l,k);
iter = 1;
%h = norm(v,2);
[Vn,Hn,vnext,hnext] = SymplecticLanczosMethod(A,v,n,conv);

invJ = [sparse(n,n),-speye(n);speye(n),sparse(n,n)];
J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
F = invJ*Vn'*J*v*F;
if PMint == 1
    Zn = int(Hn,F,ht);
elseif PMint == 2
    Zn = expintegrate(Hn,Hn\F(:,1),0:ht:ht*(k-1));
elseif PMint == 3
    Zn = real(myexpm(full(Hn),Hn\F(:,1),0:ht:ht*(k-1)));
end

ns = Vn*Zn;
U = U + ns;
diff = hnext;
if restart
    while diff > conv
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = SymplecticLanczosMethod(A,v,n,conv);
        F = invJ*Vn'*J*h*v*Zn(end,:);
        Zn = int(Hn,F,ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
    end
else
end
energy1 = max(energyBIG(A,Zn,vnext,hnext,ht,figvar,int));
energy2 = max(energySMALL(Hn,Vn,Zn,vnext,hnext,ht,figvar,int));
end