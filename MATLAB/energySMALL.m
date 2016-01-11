function energy = energySMALL(Hn,Vn,Zn,vnext,hnext,ht,figvar,int)
[n,k] = size(Zn);
l = size(Vn,1);

invJ = [sparse(n/2,n/2),-speye(n/2);speye(n/2),sparse(n/2,n/2)];
J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
Jn = [sparse(n/2,n/2),speye(n/2);-speye(n/2),sparse(n/2,n/2)];
F = invJ*Vn'*J*hnext*vnext*Zn(end,:);

delta = int(Hn,F,ht);

energy = zeros(1,k);
for i = 1:k
    energy(i) = 0.5*delta(:,i)'*Jn*Hn*delta(:,i) + delta(:,i)'*Vn'*J*hnext*vnext*Zn(end,i);
end
if figvar
    figure(3); plot(0:ht:ht*(k-1),energy,'k:.');
end
end