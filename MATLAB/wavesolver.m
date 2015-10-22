clear
close all
m = 11;
k = 100;
n = 2*(m-2);
solmeth = 2;
prob = 2;
conv = 10^-5;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k);
ht = T(2)-T(1);

%%%% Feiler er noe rare greier!!! %%%%

disk = ht^2/(hs^2);

%%% Kode starter her! %%%

%%%% Problem data %%%%

[ U0,V0,F1,F2,G1,G2,correctsolution] = getWaveTestFunctions( prob,m,k,X,T );

%%%% Matriser og vectorer %%%%
A = -1/hs^2*gallery('poisson', m-2);
Atilde = [sparse((m-2)^2,(m-2)^2),A;speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];


v = helpvector(m);
U0tilde = [U0(v);V0(v)];

F1tilde = [F1(v);sparse(length(v),1)];
F2tilde = [F2(v);sparse(length(v),1)];


if solmeth == 1
    Utemp2 = sparse(length(F1tilde),k); iter2 = 0;
    Utemp3 = sparse(length(F1tilde),k); iter3 = 0;
    tic;
    [Utemp1,iter1] = KPMwave(Atilde,Atilde*U0tilde,1,m,n,k,ht,conv,0);
    if ~(max(G1) == 0 &&  min(G1) == 0) || ~(max(F1) == 0 &&  min(F1) == 0)
    [Utemp2,iter2] = KPMwave(Atilde,F1tilde,G1,m,n,k,ht,conv,0);
    end
    if ~(max(G2) == 0 &&  min(G2) == 0) || ~(max(F2) == 0 &&  min(F2) == 0)
    [Utemp3,iter3] = KPMwave(Atilde,F2tilde,G2,m,n,k,ht,conv,0);
    end
    Utemp = Utemp1 + Utemp2 + Utemp3;
    utdata(1) = max([iter1,iter2,iter3]);
    utdata(2) = toc;
elseif solmeth == 2
    Utemp2 = sparse(length(F1tilde),k); iter2 = 0;
    Utemp3 = sparse(length(F1tilde),k); iter3 = 0;
    tic;
    [Utemp1,iter1] = KPMwave(Atilde,Atilde*U0tilde,1,m,n,k,ht,conv,1);
    if ~(max(G1) == 0 &&  min(G1) == 0) || ~(max(F1) == 0 &&  min(F1) == 0)
    [Utemp2,iter2] = KPMwave(Atilde,F1tilde,G1,m,n,k,ht,conv,1);
    end
    if ~(max(G2) == 0 &&  min(G2) == 0) || ~(max(F2) == 0 &&  min(F2) == 0)
    [Utemp3,iter3] = KPMwave(Atilde,F2tilde,G2,m,n,k,ht,conv,1);
    end
    Utemp = Utemp1 + Utemp2 + Utemp3;
    utdata(1) = max([iter1,iter2,iter3]);
    utdata(2) = toc;
elseif solmeth == 3
    utdata(1) = 1;
    tic;
    Utemp = integrate(Atilde,Atilde*U0tilde*ones(1,k)+F1tilde*G1+F2tilde*G2,2*(m-2)^2,k,ht);
    
    utdata(2) = toc;
end


U = zeros(m^2,k);
U(v,:) = Utemp(1:(m-2)^2,:);
U(v,:) = U(v,:) + U0(v)*ones(1,k);

V = zeros(m^2,k);
V(v,:) = Utemp((m-2)^2+1:end,:);
V(v,:) = V(v,:) + U0(v)*ones(1,k);

if 1
    video(U,m,k,0.05)
    video(V,m,k,0.05)
    video(correctsolution,m,k,0.05)
    video(U-correctsolution,m,k,0.05)
end

utdata(3) = max(max(max(abs(U-correctsolution))));
utdata


