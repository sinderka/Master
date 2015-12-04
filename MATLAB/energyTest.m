function utdata = energyTest(m,n,k,eqn,alg,integrator,restart,prob,conv,~)

%%% I denne koden er det masse feil!

%function utdata = energyTest
%%% Funke itj!
% Initsiell
if nargin < 10
    m = 40;
    k = 40;
    n = 4;%2*(m-2)^2-2;
    restart = 0;
    prob = 1;
    conv = 10^-14;
    para = 4; %%%%% ARG %%%%%%
    eqn = 'wave';
    alg = 2;
    integrator = 1;
end
if alg == 2
    algo = @SymplecticLanczosMethod; n = n/2;
elseif alg == 1
    algo = @Arnoldi;
end
utdata = zeros(1,6); % Burde legge til forskjellen mellom energi og
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);
[vec,height] = helpvector(m,eqn);
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V(:,1) = A*V(:,1);
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
%Utemp = 0;
[Utemp,iter1,vnext,Zn] = KPMloc(A,V(:,1),F(1,:),n,ht,conv,restart,algo,int);
%Utemp = Utemp + Utemp1;
iter = max(iter1,iter);
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
%Utemp1 = integrateloc(A,V(:,1)*F(1,:),T);
Utemp1 = expintegrate(A,U0,k,ht);
Time = toc;


%Utemp1 = Utemp1 + U0*ones(1,k);
U1 = zeros(height,k);

U1(vec,:) = Utemp1(1:length(A)/2,:);


utdata(5) = max(max(abs(U-U1)));
    %utdata(6) = abs(energy(A,Utemp-Utemp1));
    %er = Utemp-Utemp1;
utdata(6) = energy(A,Utemp-Utemp1,T,alg,Zn,vnext);
    %J = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];
    %e2n = zeros(size(Zn,1),1); e2n(end) = 1;
    %energyerror = zeros(1,k);
    %for i = 1:k
    %    energyerror(i) = 1/2*er(:,i)'*J*A*er(:,i) + er(:,i)'*J*vnext*e2n'*Zn(:,i);
    %end
    %utdata(6) = max(abs(energyerror));


%utdata(3) = getError(U,correctsolution);
utdata(3) = max(max(abs(U-correctsolution)));
utdata(4) = energy(A,Utemp);


%figure(17);plot(energyerror)
%utdata
if 0
    video(U-U1,m,k,0.05,eqn)
    pause
    video(U1,m,k,0.05,eqn)
    pause
    %V = zeros(m^2,k);
    %V(vec,:) = Utemp((m-2)^2+1:end,:);
    %V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);
    video(U,m,k,0.05,eqn)
    %video(V,m,k,0.05)
    %video(correctsolution,m,k,0.05,eqn)
    pause
    video(U-correctsolution,m,k,0.05,eqn)
    pause
    video(U1-correctsolution,m,k,0.05,eqn)
    %energy(Jtilde*Atilde,Utemp);
end
end
function [U,iter,vnext,Zn] = KPMloc(A,v,F,n,ht,conv,restart,alg,int)
%Skriv en programdefinosjon her
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
%[Zn] = integrateloc(Hn,[h*F(1,:);sparse(length(Hn)-1,k)],T);
[Zn] = int(Hn,[h*F(1,:);sparse(length(Hn)-1,k)],k,ht);

ns = Vn*Zn;
U = U + ns;
diff = hnext;
if restart
    while diff > conv
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = alg(A,v,n,conv);
        %[Zn] = integrateloc(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],T);
        [Zn] = int(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],k,ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
    end
end
vnext = hnext*vnext;
end
