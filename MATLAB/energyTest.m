function utdata = energyTest(m,n,k,eqn,alg,integrator,restart,prob,conv,~)
display('Dette programmet er for tiden ødelagt')
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
if nargin < 10
    m = 20;
    k = 20;
    n = 6;%2*(m-2)^2-2;
    restart = 0;
    prob = 1;
    conv = 10^-14;
    para = 4; %%%%% ARG %%%%%%
    eqn = 'semirandom';
    alg = 2;
    integrator = 1;
end
if alg == 2
    algo = @SymplecticLanczosMethod; n = n/2;
elseif alg == 1
    algo = @Arnoldi;
end
%%% Initsiell data
utdata = zeros(1,6);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);
[vec,height] = helpvector(m,eqn);

% Get problem information
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V(:,1) = A*V(:,1);

% Chose integration method
if integrator == 1
    int = @trapezoidal;
elseif integrator == 2
    int = @forwardeuler;
elseif integrator == 3
    int = @implicitmidpoint;
elseif integrator == 4
    int = @expintegrate;
end

tic;
iter = 0;
Utemp = 0;
for i = 1:size(F,1)
    [Utemp1,iter1,Hn,Vn,Zn] = KPMloc(A,V(:,i),F(i,:),n,ht,conv,restart,algo,int);
    Utemp = Utemp + Utemp1;
    iter = max(iter1,iter);
end% End forloop


utdata(1) = iter;
utdata(2) = toc;
U = zeros(height,k);
if integrator ~= 4
    Utemp = Utemp + U0*ones(1,k);
end
U(vec,:) = Utemp(1:length(A)/2,:);

tic;
tempVF = 0;
for i = 1:size(F,1)
    tempVF = tempVF + V(:,i)*F(i,:);
end
Utemp1 = int(A,tempVF,ht);
%Utemp1 = expintegrate(A,U0,k,ht);
%Time = toc;


Utemp1 = Utemp1 + U0*ones(1,k);
U1 = zeros(height,k);

U1(vec,:) = Utemp1(1:length(A)/2,:);


utdata(5) = max(max(abs(U-U1)));
figure(5);plot(T,max(abs(U-U1)), 'k:.');

figure(7); plot(T,energy(A,Utemp,U0,alg,Hn,Vn,Zn,restart) - energy(A,Utemp1,U0),'k:.');
utdata(6) = max(abs(energy(A,Utemp,U0,alg,Hn,Vn,Zn,restart) - energy(A,Utemp1,U0)));

figure(2);plot(T,energy(A,Utemp,U0,alg,Hn,Vn,Zn,restart),'k:.')
utdata(4) = max(abs(energy(A,Utemp,U0,alg,Hn,Vn,Zn,restart)));

figure(11);plot(T,max(abs(U-correctsolution)),'k:.');
utdata(3) = max(max(abs(U-correctsolution)));


if 0
    video(U-U1,m,k,0.05,eqn)
    %pause
    video(U1,m,k,0.05,eqn)
    %pause
    %V = zeros(m^2,k);
    %V(vec,:) = Utemp((m-2)^2+1:end,:);
    %V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);
    %video(U,m,k,0.05,eqn)
    %video(V,m,k,0.05)
    %video(correctsolution,m,k,0.05,eqn)
    %pause
    %video(U-correctsolution,m,k,0.05,eqn)
    %pause
    %video(U1-correctsolution,m,k,0.05,eqn)
    %energy(Jtilde*Atilde,Utemp);
end
end
function [U,iter,Hn,Vn,Zn] = KPMloc(A,v,F,n,ht,conv,restart,alg,int)
l = size(A,1);
k = length(F);
if max(abs(v)) == 0 || max(max(abs(F))) == 0
    U = sparse(l,k);
    iter = 0;
    return
end
U = zeros(l,k);
iter = 1;
h = norm(v,2);
[Vn,Hn,vnext,hnext] = alg(A,v,n,conv);
[Zn] = int(Hn,[h*F(1,:);sparse(length(Hn)-1,k)],ht);

ns = Vn*Zn;
U = U + ns;
diff = hnext;
if restart
    while diff > conv
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = alg(A,v,n,conv);
        [Zn] = int(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
    end
end
vnext = hnext*vnext;
end
