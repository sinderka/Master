function Max_Energy_Difference = energy(A,y,T,alg,Zn,vnext)
if nargin == 3
    alg = 1;
end
%Skriv en programdefinosjon her



if alg == 1 || alg == 3
    k = size(y,2);
    m = length(A)/2;
    A = [sparse(m,m),speye(m);-speye(m),sparse(m,m)]*A;
    ener = zeros(1,k);
    initEnergy = 0.5*y(:,1)'*A*y(:,1);
    for i = 1:k
        ener(i) = 0.5*y(:,i)'*A*y(:,i);
    end
    
    figure(2); plot(T,initEnergy-ener, 'k:.')
    Max_Energy_Difference = max(abs(ener-initEnergy));
elseif alg == 2
    k = size(y,2);
    m = length(A)/2;
    J = [sparse(m,m),speye(m);-speye(m),sparse(m,m)];
    e2n = zeros(size(Zn,1),1); e2n(end) = 1;
    energyerror = zeros(1,k);
    for i = 1:k
        energyerror(i) = 1/2*y(:,i)'*J*A*y(:,i) + y(:,i)'*J*vnext*e2n'*Zn(:,i);
    end
    figure(2); plot(T,energyerror(1)-energyerror, 'k:.')
    Max_Energy_Difference = max(abs(energyerror(1)-energyerror));
    
end

