function utdata = energyTest(m,n,k,simtime,eqn,restart,prob,conv)
% Solves a problem dependant on the indata
%input
% m: number of points in eqch spacial direction X
% n: restart variable (size of orthogonal space) X
% k: number of space in time. X
% eqn: says something about with algorithm to solve 
% alg(1,2,3): declares with ortogonalisation method to use X
% integrator(1,2,3): declares with integration method to use X
% restart(0,1): should the method restart or not X
% prob: says smoething about with particular problem to solve
% conv: convergence criterion used in arnoldi og KPM X
% para: currently nothing
%returns:
% utdata:
% utdata(1): number of iteration performed
% utdata(2): computation time
% utdata(3): error
% utdata(4): energy
% utdata(5): error difference
% utdata(6): energy difference
% Initsiell
if nargin < 13
    m = 25;
    simtime = 1;
    %K = 2;
    k = 25;
    n = 4;%2*(m-2)^2;
    restart = 1;
    prob = 1;
    conv = 10^-14;
    %para = 4; %%%%% If need be %%%%%%
    eqn = 'wave';
    %alg = 2;
    %integrator = 3;
    figvar = 1;
end

utdata = zeros(1,6);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,simtime,k); ht = T(2)-T(1);
[vec,height,lastrelevant] = helpvector(m,eqn);

% Get problem information
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V(:,1) = A*V(:,1);

% Chose solution method and solve
U = zeros(2*m^2,k);

for i = 1:size(F,1)
    [Utemp] = KPMloc(A,V(:,i),T,n/2,conv,restart,ht,figvar);
    U(vec,:) = U(vec,:) -Utemp;
end
U(vec,:) = U(vec,:) + U0*ones(1,k);




if 1
    %video(U(1:m^2,:),m,0.05,eqn)
    %video(correctsolution(1:m^2,:),m,0.05,eqn)
    video(U(1:m^2,:)-correctsolution(1:m^2,:),m,0.05,eqn)
end
end

function [U] = KPMloc(A,v,T,n,conv,restart,ht,figvar)
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
k = length(T);
if max(abs(v)) == 0
    U = sparse(l,k);
    Zn = sparse(2*n,k);
    vnext = sparse(l,1);
    hnext = sparse(1,1);
    %iter = 0;
    return
end
U = zeros(l,k);
%iter = 1;
h = norm(v,2);
[Vn,Hn,vnext,hnext] = SymplecticLanczosMethod(A,v,n,conv);

[Zn] = intloc(Hn,Vn,v,T);

ns = Vn*Zn;
U = U + ns;
diff = hnext;
if restart
    
    invJ = [sparse(n,n),-speye(n);speye(n),sparse(n,n)];
    J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
    
    while diff > conv 
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = SymplecticLanczosMethod(A,v,n,conv);
        %[Zn] = trapezoidal(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],ht); %Hvordan blir restarten med SLM??
        F = invJ*Vn'*J*h*v*Zn(end,:);
        Zn = trapezoidal(Hn,F,ht);
        %[Zn] = intloc(Hn,Vn,v,T);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        %iter = iter+1;
    end
else
energyBIG(A,Zn,vnext,hnext,ht,figvar)
energySMALL(Hn,Vn,Zn,vnext,hnext,ht,figvar)
end

end

function z = intloc(H,S,b,T)
[m,n] = size(S);
z = zeros(length(H),length(T));
invJ = [sparse(n/2,n/2),-speye(n/2);speye(n/2),sparse(n/2,n/2)];
J = [sparse(m/2,m/2),speye(m/2);-speye(m/2),sparse(m/2,m/2)];
invJSJb = invJ*S'*J*b;
for i = 1:length(T)
    z(:,i) = (speye(n)-expm(H*T(i)))*(H\invJSJb);
end
a = 2;
end

function energyBIG(A,Zn,vnext,hnext,ht,figvar)
m = length(A); [n,k] = size(Zn);
F = hnext*vnext*Zn(end,:);


%en = zeros(n,1); en(end) = 1;
%hnext*vnext*en'*Zn


epsilon = trapezoidal(A,F,ht);

Jm = [sparse(m/2,m/2),speye(m/2);-speye(m/2),sparse(m/2,m/2)];
energyBIG = zeros(1,k);
for i = 1:k
    energyBIG(i) = 0.5*epsilon(:,i)'*Jm*A*epsilon(:,i)+epsilon(:,i)'*Jm*hnext*vnext*Zn(end,i);
end
if figvar
    figure(2); plot(0:ht:ht*(k-1),energyBIG,'k:.');
end
end

function energySMALL(Hn,Vn,Zn,vnext,hnext,ht,figvar)
[n,k] = size(Zn);
l = size(Vn,1);

invJ = [sparse(n/2,n/2),-speye(n/2);speye(n/2),sparse(n/2,n/2)];
J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
Jn = [sparse(n/2,n/2),speye(n/2);-speye(n/2),sparse(n/2,n/2)];
F = invJ*Vn'*J*hnext*vnext*Zn(end,:);

delta = trapezoidal(Hn,F,ht);

energySMALL = zeros(1,k);
for i = 1:k
    energySMALL(i) = 0.5*delta(:,i)'*Jn*Hn*delta(:,i) + delta(:,i)'*Vn'*J*hnext*vnext*Zn(end,i);
end
if figvar
    figure(3); plot(0:ht:ht*(k-1),energySMALL,'k:.');
end
end