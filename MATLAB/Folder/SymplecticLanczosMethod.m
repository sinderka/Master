function [S,Htilde,Vend,xiend] = SymplecticLanczosMethod(H,J,v,var)

n = length(H)/2;

delta = zeros(var,1);
beta = zeros(var,1);
xi = zeros(var+1,1);
nu = zeros(var,1);

V = zeros(2*n,var+1);
W = zeros(2*n,var);


xi(2) = norm(v,2);

V(:,2) = 1/xi(2)*v;
%%%% PROBLEMER %%%%
% S'*J*S ~= J,  S'*J*S == -J = J'


for m = 1:1:var
   % Computing v
   v = H*V(:,m+1);
   % Computing delta
   delta(m) = V(:,m+1)'*v;
   % Computing Wm
   wtilde = v-delta(m)*V(:,m+1);
   
   
   %nu(m) = v'*J*V(:,m+1);
   %nu(m) = v'*inv(J)*V(:,m+1);
   nu(m) = V(:,m+1)'*J*v; %Fungerer ish!
   %nu(m) = V(:,m+1)'*inv(J)*v;
   
   
   W(:,m) = 1/nu(m)*wtilde;
   % Computing w
   w = H*W(:,m);
   % Computing beta
   
   %beta(m) = -w'*J*W(:,m);
   %beta(m) = -w'*inv(J)*W(:,m);
   beta(m) = -W(:,m)'*J*w; %Fungerer ish!
   %beta(m) = -W(:,m)'*inv(J)*w;
   
   
   
   %Computing Wm+1
   vmtilde = w-xi(m+1)*V(:,m)-beta(m)*V(:,m+1)+delta(m)*W(:,m);
   xi(m+2) = norm(vmtilde,2);
   V(:,m+2) = 1/xi(m+2)*vmtilde;
end
S = [V(:,2:end-1),W];

Htilde = [sparse(1:var,1:var,delta,var,var),gallery('tridiag',xi(3:end-1),beta,xi(3:end-1));
          sparse(1:var,1:var,nu,var,var), sparse(1:var,1:var,-delta,var,var)];


Vend = V(:,end); xiend = xi(end);
end

