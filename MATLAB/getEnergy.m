function energychange = getEnergy(A,y,b)
% returns the energy to the system
%Indata
% A: an mxm matrix
% y: a matrix of the solution
% T: a list of points in time
%Optional inputs(only for slm without restart):
% alg: numerical values of 1,2,3, depending on the orthogonalisation method used.
% Zn: return from slm
% vnext: residual vector


if nargin == 2
    b = zeros(length(A),1);
    
end


k = size(y,2);
m = length(A)/2;
Jm = [sparse(m,m),speye(m);-speye(m),sparse(m,m)];
energyerror = zeros(1,k);
for i = 1:k
    energyerror(i) = 0.5*y(:,i)'*Jm*A*y(:,i)+ y(:,i)'*Jm*b;
end
energychange = abs(energyerror(1)-energyerror);

end

