clear
close all
%%% Initsiell data
m = 40;
k = 40;
n = 40;
solmeth = 2;
prob = 5;
conv = 10^-5;
para = 4; %%%%% ARG %%%%%%
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);

vec = helpvector(m);

disk = ht^2/(hs^2);

%%%%% TODO %%%%%
%%%% Fjerne if løkkene fra programmet.
%%%% Skrive Restart og Krylov sammen


[ U0,V0,F1,F2,G1,G2,correctsolution] = getWaveTestFunctions( prob,m,k,X,T );

%%%% Matriser og vectorer %%%%
A = -1/hs^2*gallery('poisson', m-2);
%Atilde = [sparse((m-2)^2,(m-2)^2),A;speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];
Atilde = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);A,sparse((m-2)^2,(m-2)^2)];
U0tilde = [U0(vec);V0(vec)];
%F1tilde = [F1(vec);sparse((m-2)^2,1)];
%F2tilde = [F2(vec);sparse((m-2)^2,1)];
F1tilde = [sparse((m-2)^2,1);F1(vec)];
F2tilde = [sparse((m-2)^2,1);F2(vec)];


% A = -1/hs^2*gallery('poisson', m); %%% FIX???
% %Atilde = [sparse((m-2)^2,(m-2)^2),A;speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];
% Atilde = [sparse(m^2,m^2),speye(m^2);A,sparse(m^2,m^2)];
% U0tilde = [U0;V0];
% F1tilde = [F1;sparse(m^2,1)];
% F2tilde = [F2;sparse(m^2,1)];




if solmeth == 1
    Utemp1 = sparse(2*(m-2)^2,k); iter1 = 0;
    Utemp2 = sparse(2*(m-2)^2,k); iter2 = 0;
    Utemp3 = sparse(2*(m-2)^2,k); iter3 = 0;
    tic;
    if ~(max(U0tilde) == 0 &&  min(U0tilde) == 0)
        [Utemp1,iter1] = KPM(Atilde,Atilde*U0tilde,1,k,n,2*(m-2)^2,ht,conv,0);
    end
    if ~(max(G1) == 0 &&  min(G1) == 0) || ~(max(F1) == 0 &&  min(F1) == 0)
        [Utemp2,iter2] = KPM(Atilde,F1tilde,G1,k,n,2*(m-2)^2,ht,conv,0);
    end
    if ~(max(G2) == 0 &&  min(G2) == 0) || ~(max(F2) == 0 &&  min(F2) == 0)
        [Utemp3,iter3] = KPM(Atilde,F2tilde,G2,k,n,2*(m-2)^2,ht,conv,0);
    end
    Utemp = Utemp1 + Utemp2 + Utemp3;
    utdata(1) = max([iter1,iter2,iter3]);
    utdata(2) = toc;
elseif solmeth == 2
    Utemp1 = sparse(length(F1tilde),k); iter1 = 0;
    Utemp2 = sparse(length(F1tilde),k); iter2 = 0;
    Utemp3 = sparse(length(F1tilde),k); iter3 = 0;
    tic;
    if ~(max(U0tilde) == 0 &&  min(U0tilde) == 0)
        [Utemp1,iter1] = KPM(Atilde,Atilde*U0tilde,1,k,n,2*(m-2)^2,ht,conv,1);
    end
    if ~(max(G1) == 0 &&  min(G1) == 0) || ~(max(F1) == 0 &&  min(F1) == 0)
        [Utemp2,iter2] = KPM(Atilde,F1tilde,G1,k,n,2*(m-2)^2,ht,conv,1);
    end
    if ~(max(G2) == 0 &&  min(G2) == 0) || ~(max(F2) == 0 &&  min(F2) == 0)
        [Utemp3,iter3] = KPM(Atilde,F2tilde,G2,k,n,2*(m-2)^2,ht,conv,1);
    end
    Utemp = Utemp1 + Utemp2 + Utemp3;
    utdata(1) = max([iter1,iter2,iter3]);
    utdata(2) = toc;
elseif solmeth == 3
    utdata(1) = 0;
    tic;
    Utemp = integrate(Atilde,Atilde*U0tilde*ones(1,k)+F1tilde*G1+F2tilde*G2,2*(m-2)^2,k,ht);
    utdata(2) = toc;
end

U = zeros(m^2,k);
U(vec,:) = Utemp(1:(m-2)^2,:);
U(vec,:) = U(vec,:) + U0(vec)*ones(1,k);

V = zeros(m^2,k);
V(vec,:) = Utemp((m-2)^2+1:end,:);
V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);

if 0
    video(U,m,k,0.05)
    %video(V,m,k,0.05)
    video(correctsolution,m,k,0.05)
    video(U-correctsolution,m,k,0.05)
end

utdata(3) = max(max(abs(U-correctsolution)));
utdata
%utdata(3)

