clear
close all
%%% Initsiell data
m = 20;
k = 20;
n = (m-2)^2;
prob = 1;
solmeth = 3;
conv = 10^-15;
para = 4; %%%%% ARG %%%%%%
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m); hs = X(2) - X(1);
T = linspace(0,1,k); ht = T(2) - T(1);

%%%%%%%%% TODO %%%%%%%%%
%%%% Skrive Restart og Krylov sammen
%%%% Finne ut hvorfor Restart og Krylov ikke konvergerer
%%%% Fjerne if løkkene fra programmet.

vec = helpvector(m);

disk = ht^2/(hs^2);

[U0,F1,F2,G1,G2,correctsolution] = getHeatTestFunction(prob,m,k,X,T);

%%%% Matriser og vectorer %%%%
A = -1/hs^2*gallery('poisson', m-2);
%A = -1/hs^2*gallery('poisson', m); %%% FIX???

U = zeros(m^2,k);
U1 = sparse((m-2)^2,k); iter1 = 0;
U2 = sparse((m-2)^2,k); iter2 = 0;
U3 = sparse((m-2)^2,k); iter3 = 0;
if solmeth == 1 %% Krylov
    tic;
    %[U1,iter1] = KPM(A,A*U0(vec),1,k,n,(m-2)^2,ht,conv,0);
    [U2,iter2] = KPM(A,F1(vec),G1,k,n,(m-2)^2,ht,conv,0);
    [U3,iter3] = KPM(A,F2(vec),G2,k,n,(m-2)^2,ht,conv,0);
    U(vec,:) = U1 + U2 + U3;
    utdata(1) = max([iter1,iter2,iter3]);
    utdata(2) = toc;
elseif solmeth == 2 %% Restarted krylov
    tic;
    %[U1,iter1] = KPM(A,A*U0(vec),1,k,n,(m-2)^2,ht,conv,1);
    [U2,iter2] = KPM(A,F1(vec),G1,k,n,(m-2)^2,ht,conv,1);
    [U3,iter3] = KPM(A,F2(vec),G2,k,n,(m-2)^2,ht,conv,1);
    U(vec,:) = U1 + U2 + U3;
    utdata(1) = max([iter1,iter2,iter3]);
    utdata(2) = toc;
elseif solmeth == 3 %% Direkte integrasjon
    utdata(1) = 0;
    tic;
    U(vec,:) = integrate(A,A*U0(vec)*ones(1,k)+F1(vec)*G1+F2(vec)*G2,(m-2)^2,k,ht);
    %U = integrate(A,F1*G1+F2*G2,m^2,k,ht);
    utdata(2) = toc;
end
U(vec,:) = U(vec,:) + U0(vec)*ones(1,k);

if 0
    video(U,m,k,0.05)
    video(correctsolution,m,k,0.05)
    video(U-correctsolution,m,k,0.05)
end

utdata(3) = max(max(abs(U-correctsolution)));
utdata
%utdata(3)