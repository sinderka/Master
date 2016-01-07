function energychange = energy(A,y,alg,Hn,Vn,Zn,restart)
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
    alg = 1;
end



if 1%alg == 1 || alg == 3 || restart == 1
    k = size(y,2);
    m = length(A)/2;
%     if size(y,1) ~= 2*m
%         m = (m-1)/2;
%     end
    %A = [sparse(m,m),speye(m);-speye(m),sparse(m,m)]*A;
    Jm = [sparse(m,m),speye(m);-speye(m),sparse(m,m)];
    energyerror = zeros(1,k);
    for i = 1:k
        %energyerror(i) = 0.5*y(:,i)'*Jm*A*y(:,i) + y(:,i)'*Jm*U0;
        energyerror(i) = 0.5*y(:,i)'*Jm*A*y(:,i);%+ y(:,i)'*Jm*U0;% + y(:,i)'*Jm*A*U0;
    end
    energychange = energyerror(1)-energyerror;
elseif alg == 2 && restart == 0
    k = size(y,2);
    m = length(A)/2;
    n = size(Hn,2)/2;
    Jm = [sparse(m,m),speye(m);-speye(m),sparse(m,m)];
    Jn = [sparse(n,n),speye(n);-speye(n),sparse(n,n)];
    energyerror = zeros(1,k);
    for i = 1:k
        energyerror(i) = 0.5*Zn(:,i)'*Jn*Hn*Zn(:,i) + Zn(:,i)'*Vn'*Jm*U0;
        %energyerror(i) = 1/2*y(:,i)'*J*A*y(:,i) + y(:,i)'*J*vnext*e2n'*Zn(:,i);
    end
    energychange = energyerror(1)-energyerror;
    
end

