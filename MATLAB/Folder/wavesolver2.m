clear
close all
m = 11;
k = 100;
n = 1;
solmeth = 2;
prob = 1;
conv = 10^-5;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k);
ht = T(2)-T(1);
A = -1/hs^2*gallery('poisson', m-2);

disk = ht^2/(hs^2);

%%% Kode starter her! %%%

%%%% Problem data %%%%

[ U0,V0,F,correctsolution] = getWaveTestFunctions( prob,m,k,X,T );

%%%% Matriser og vectorer %%%%

Atilde = [sparse((m-2)^2,(m-2)^2),A;speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];
v = helpvector(m);
U0tilde = [U0(v);V0(v)];
%U0tilde = [V0(v);U0(v)];

U = zeros(m^2,k);
%V = zeros(m^2,k);
if solmeth == 1
    tic;
    [Utemp,iter] = KPMwave2(Atilde,U0tilde,m,2*(m-2)^2,k,ht,conv);
    utdata(1) = iter;
    utdata(2) = toc;
elseif solmeth == 2
    tic;
    [Utemp,iter] = KPMwave2(Atilde,U0tilde,m,n,k,ht,conv);
    utdata(1) = iter;
    utdata(2) = toc;
elseif solmeth == 3
    utdata(1) = 1;
    tic;
    Utemp = integrate2(Atilde,Atilde*U0tilde,2*(m-2)^2,k,ht);
    
    utdata(2) = toc;
end
U(v,:) = Utemp(1:(m-2)^2,:);
U(v,:) = U(v,:) + U0(v)*ones(1,k);
if 1
    video(U,m,k,0.05)
    %video(V,m,k,0.05)
    video(correctsolution,m,k,0.05)
    video(U-correctsolution,m,k,0.05)
end

utdata(3) = max(max(max(abs(U-correctsolution))));
utdata


