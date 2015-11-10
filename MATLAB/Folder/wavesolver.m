%Skriv en programdefinosjon her
clear
close all
%%% Initsiell data
%tic;
m = 30;
k = 30;
n = 4;%2*(m-2)^2;
solmeth = 1;
restart = 1;
prob = 7;
conv = 10^-10;
para = 4; %%%%% ARG %%%%%%
eqn = 'wave';
alg = 1;
%Alt over dette burde være argumenter +

%%% Initsiell data
utdata = zeros(1,4);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);

vec = helpvector(m);

disk = ht^2/(hs^2);


%%%%%%%%%% TODO %%%%%%%%%%
%%% Skrive bølge, maxwell og varme sammen så mye som mulig %%%
%%% Legge til beskrivelse til ferdige funksjoner %%%
%%%


[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getWaveTestFunctions( prob,X,T ); %Denne skal byttes ut med "getTestFunctions(prob,X,T,var)"
V(:,1) = A*V(:,1); % denne kodelinjen er teit?
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAN:
% funksjon
% Initial stuff X
% Verifiser inndata
% Få tak i testfunksjoner X--
% en forløkke som løser problemet X(-restart for SLM)
% beregne feil
% annet

if alg == 1
    alg = @Arnoldi;
elseif alg == 2
    alg = @SymplecticLanczosMethod;
end
if solmeth == 1
    tic;
    iter = 0;
    Utemp = 0;
    for i = 1:size(F,1)
        [Utemp1,iter1] = KPM(A,V(:,i),F(i,:),n,ht,conv,restart,alg);
        Utemp = Utemp + Utemp1;
        iter = max(iter1,iter);
    end
    utdata(1) = iter;
    utdata(2) = toc;
elseif solmeth == 2
    utdata(1) = 0;
    tic;
    tempVF = 0;
    for i = 1:size(F,1)
        tempVF = tempVF + V(:,i)*F(i,:);
    end
    Utemp = integrate(A,tempVF,k,ht);
    utdata(2) = toc;
    
end
U = zeros(m^2,k);
Utemp = Utemp + U0*ones(1,k);
U(vec,:) = Utemp(1:(m-2)^2,:);

utdata(3) = getError(U,correctsolution);

if 1
    %V = zeros(m^2,k);
    %V(vec,:) = Utemp((m-2)^2+1:end,:);
    %V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);
    video(U,m,k,0.05)
    %video(V,m,k,0.05)
    video(correctsolution,m,k,0.05)
    video(U-correctsolution,m,k,0.05)
    %energy(Jtilde*Atilde,Utemp);
end
utdata(4) = energy(A,Utemp);
%energy(Jtilde*Atilde,Utemp-correctsolution(vec,(m)^2+vec,:)); % Fungerer
%ikke fordi den deriverte ikke blir tatt med i beregningene av energi


error = max(max(abs(U-correctsolution)))
utdata
%utdata(3)
%toc;