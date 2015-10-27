clear
close all
%%% Initsiell data
m = 41;
k = 41;
n = 41;
solmeth = 1;
restart = 1;
prob = 5;
conv = 10^-5;
para = 4; %%%%% ARG %%%%%%
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);

vec = helpvector(m);

disk = ht^2/(hs^2);

[ U0,V0,F1,F2,G1,G2,correctsolution] = getWaveTestFunctions( prob,m,k,X,T );

%%%% Matriser og vectorer %%%%
A = 1/hs^2*gallery('poisson', m-2);
Atilde = [A,sparse((m-2)^2,(m-2)^2);sparse((m-2)^2,(m-2)^2),speye((m-2)^2)]; %% Korrekt
Jtilde = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];
    
%Atilde = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);A,sparse((m-2)^2,(m-2)^2)]; %Test
%Jtilde = [speye((m-2)^2),sparse((m-2)^2,(m-2)^2);sparse((m-2)^2,(m-2)^2),-speye((m-2)^2)];

U0tilde = [U0(vec);V0(vec)];
F1tilde = [sparse((m-2)^2,1);F1(vec)];
F2tilde = [sparse((m-2)^2,1);F2(vec)];

if solmeth == 1
    tic;
    [Utemp1,iter1] = KPM(Jtilde,Atilde,Jtilde*Atilde*U0tilde,1,k,n,2*(m-2)^2,ht,conv,restart);
    [Utemp2,iter2] = KPM(Jtilde,Atilde,F1tilde,G1,k,n,2*(m-2)^2,ht,conv,restart);
    [Utemp3,iter3] = KPM(Jtilde,Atilde,F2tilde,G2,k,n,2*(m-2)^2,ht,conv,restart);
    Utemp = Utemp1 + Utemp2 + Utemp3;
    utdata(1) = iter1+iter2+iter3;
    utdata(2) = toc;
elseif solmeth == 2
    utdata(1) = 0;
    tic;
    Utemp = integrate(Jtilde*Atilde,Jtilde*Atilde*U0tilde*ones(1,k)+F1tilde*G1+F2tilde*G2,2*(m-2)^2,k,ht);
    utdata(2) = toc;
end

U = zeros(m^2,k);
U(vec,:) = Utemp(1:(m-2)^2,:);
U(vec,:) = U(vec,:) + U0(vec)*ones(1,k);

%V = zeros(m^2,k);
%V(vec,:) = Utemp((m-2)^2+1:end,:);
%V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);

if 1
    video(U,m,k,0.05)
    %video(V,m,k,0.05)
    video(correctsolution,m,k,0.05)
    video(U-correctsolution,m,k,0.05)
end

utdata(3) = max(max(abs(U-correctsolution)));
utdata
%utdata(3)