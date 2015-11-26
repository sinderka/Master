function [utdata,energyerror] =  energyTest
%%% Funke itj!
% Initsiell
m = 20;
k = 20;
n = 8;%2*(m-2)^2;
restart = 0;
prob = 1;
conv = 10^-14;
para = 4; %%%%% ARG %%%%%%
eqn = 'semirandom';
alg = 1;
utdata = zeros(1,6); % Burde legge til forskjellen mellom energi og
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);
[vec,height] = helpvector(m,eqn);
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V(:,1) = A*V(:,1);
if alg == 1 || alg == 2
    if alg == 1
        algo = @Arnoldi;
    elseif alg == 2
        algo = @SymplecticLanczosMethod;
        n = n/2;
    end
    tic;
    iter = 0;
    Utemp = 0;
    for i = 1:1%size(F,1)
        [Utemp1,iter1,vnext,Zn] = KPMloc(A,V(:,i),F(i,:),n,ht,conv,restart,algo,T);
        Utemp = Utemp + Utemp1;
        iter = max(iter1,iter);
    end
    utdata(1) = iter;
    utdata(2) = toc;
    U = zeros(height,k);
    Utemp = Utemp + U0*ones(1,k);
    U(vec,:) = Utemp(1:length(A)/2,:);
end
tic;
tempVF = 0;
for i = 1:size(F,1)
    tempVF = tempVF + V(:,i)*F(i,:);
end
%Utemp1 = integrateloc(A,V(:,1)*F(1,:),T);
Utemp1 = integrate(A,V(:,1)*F(1,:),k,ht);
Time = toc;


Utemp1 = Utemp1 + U0*ones(1,k);
U1 = zeros(height,k);
if strcmp(eqn,'maxwell1D')
    U1(vec,:) = Utemp1(1:length(A)/2-1,:); % OBS: Dette er en dårlig løsning!
else
    U1(vec,:) = Utemp1(1:length(A)/2,:);
end


if alg ~= 3
    utdata(5) = max(max(abs(U-U1)));
    utdata(6) = abs(energy(A,Utemp-Utemp1));
else
    utdata(1) = 0;
    utdata(2) = Time;
    utdata(5) = -1;
    utdata(6) = -1;
    Utemp = Utemp1;
    U = U1;
end


%utdata(3) = getError(U,correctsolution);
utdata(3) = max(max(abs(U-correctsolution)));
utdata(4) = energy(A,Utemp);

er = Utemp-Utemp1;
J = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];
e2m = zeros(size(Zn,1),1); e2m(end) = 1;
energyerror = zeros(1,k);
for i = 1:k
    energyerror(i) = 1/2*er(:,i)'*J*A*er(:,i) + er(:,i)'*J*vnext*e2m'*Zn(:,i);
end
utdata
if 0
    %V = zeros(m^2,k);
    %V(vec,:) = Utemp((m-2)^2+1:end,:);
    %V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);
    video(U,m,k,0.05,eqn)
    %video(V,m,k,0.05)
    %video(correctsolution,m,k,0.05,eqn)
    %video(U-correctsolution,m,k,0.05,eqn)
    %energy(Jtilde*Atilde,Utemp);
end
end
function [U,iter,vnext,Zn] = KPMloc(A,v,F,n,ht,conv,restart,alg,T)
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
[Zn] = integrate(Hn,[h*F(1,:);sparse(length(Hn)-1,k)],k,ht);

ns = Vn*Zn;
U = U + ns;
diff = hnext;
if restart
    while diff > conv
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = alg(A,v,n,conv);
        %[Zn] = integrateloc(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],T);
        [Zn] = integrate(Hn,[h*Zn(end,:);sparse(length(Hn)-1,k)],k,ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
    end
end
vnext = hnext*vnext;
end
function U = integrateloc(A,U0,T)
if size(U0,2) == 1
    U = zeros(length(U0),length(T));
    for i = 1:length(T)
        U(:,i) = exp(A*T(i))*U0;
    end
else
    U = zeros(length(U0),length(T));
    for i = 1:length(T)
        U(:,i) = exp(A*T(i))*U0(:,i);
    end
end
end
