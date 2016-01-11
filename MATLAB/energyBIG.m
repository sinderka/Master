function energy = energyBIG(A,Zn,vnext,hnext,ht,figvar,int)
m = length(A); [n,k] = size(Zn);
F = hnext*vnext*Zn(end,:);


%en = zeros(n,1); en(end) = 1;
%hnext*vnext*en'*Zn


epsilon = int(A,F,ht);

Jm = [sparse(m/2,m/2),speye(m/2);-speye(m/2),sparse(m/2,m/2)];
energy = zeros(1,k);
for i = 1:k
    energy(i) = 0.5*epsilon(:,i)'*Jm*A*epsilon(:,i)+epsilon(:,i)'*Jm*hnext*vnext*Zn(end,i);
end
if figvar
    figure(2); plot(0:ht:ht*(k-1),energy,'k:.');
end

end