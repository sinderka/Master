function Max_Energy_Difference = energy(A,y)
%Skriv en programdefinosjon her
k = size(y,2);
m = length(A)/2;
A = [sparse(m,m),speye(m);-speye(m),sparse(m,m)]*A;

ener = zeros(1,k);
initEnergy = 0.5*y(:,1)'*A*y(:,1);
for i = 1:k
ener(i) = 0.5*y(:,i)'*A*y(:,i);
end

figure(2)
plot(1:k,initEnergy-ener)
Max_Energy_Difference = max(max(ener))-min(min(ener));

end

