function [S,Htilde,Vend,xiend] = SymplecticLanczosMethod(H,J,v,var)
%Skriv en programdefinosjon her
n = length(H)/2;

delta = zeros(var,1);
beta = zeros(var,1);
xi = zeros(var+1,1);
nu = zeros(var,1);

V = zeros(2*n,var+1);
W = zeros(2*n,var);


xi(2) = norm(v,2);

V(:,2) = 1/xi(2)*v;

for m = 1:1:var
   % Computing v
   v = H*V(:,m+1);
   % Computing delta
   delta(m) = V(:,m+1)'*v;
   % Computing Wm
   wtilde = v-delta(m)*V(:,m+1);
   
   nu(m) = V(:,m+1)'*J*v; 
   
   W(:,m) = 1/nu(m)*wtilde;
   W(:,m) = W(:,m)+[V(:,2:m),W(:,1:m-1)]*[sparse(m-1,m-1),speye(m-1);-speye(m-1),sparse(m-1,m-1)]*[V(:,2:m),W(:,1:m-1)]'*J*W(:,m);
   % Computing w
   w = H*W(:,m);
 
   % Computing beta
   beta(m) = -W(:,m)'*J*w;
   %Computing Wm+1
   vmtilde = w-xi(m+1)*V(:,m)-beta(m)*V(:,m+1)+delta(m)*W(:,m);
   xi(m+2) = norm(vmtilde,2);
   V(:,m+2) = 1/xi(m+2)*vmtilde;
   V(:,m+2) = V(:,m+2)+[V(:,2:m+1),W(:,1:m)]*[sparse(m,m),speye(m);-speye(m),sparse(m,m)]*[V(:,2:m+1),W(:,1:m)]'*J*V(:,m+2);
end

S = [V(:,2:end-1),W];

Htilde = [sparse(1:var,1:var,delta,var,var),gallery('tridiag',xi(3:end-1),beta,xi(3:end-1));
          sparse(1:var,1:var,nu,var,var), sparse(1:var,1:var,-delta,var,var)];


Vend = V(:,end); xiend = xi(end);
end

