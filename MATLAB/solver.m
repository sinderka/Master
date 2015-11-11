function utdata = solver(m,n,k,eqn,alg,restart,prob,conv,para)
%Rekkefølge inn: m,n,k,eqn,alg,restart,prob,conv,para
%Skriv en programdefinosjon her
%clear
%close all
%%% Initiell data
%tic;
if nargin < 9
    m = 40;
    k = 10;
    n = 4;%2*(m-2)^2;
    solmeth = 1;
    restart = 1;
    prob = 1;
    conv = 10^-14;
    para = 4; %%%%% ARG %%%%%%
    eqn = 'wave';
    alg = 1;
end
%Alt over dette burde være argumenter +
solmeth = 1;
if ~verifyData(m,n,k,eqn,alg,restart,prob,para)
    utdata = -ones(1,4);
    return
end
%m,n,k,eqn,alg,restart,prob,conv,para
%datastore1 = datastorage(m,n,k,restart,eqn,prob,para,alg,conv,1) % Sjekker om krylov må
%benyttes
%datastore2 = datastorage(m,n,k,restart,eqn,prob,para,alg,conv,2) % Sjekker
%om den er løst av direkte integrasjons metoden


%%% Initsiell data
utdata = zeros(1,5); % Burde legge til forskjellen mellom energi og
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);

vec = helpvector(m);

%disk = ht^2/(hs^2);


%%%%%%%%%% TODO %%%%%%%%%%
%%% Skrive bølge, maxwell og varme sammen så mye som mulig %%%
%%% Legge til beskrivelse til ferdige funksjoner %%%
%%% Lage en database av resultater, så beregninger går fortere!
%%% Feilen mellom direkte metode og Krylov metode burde også være med,
%%% alltid!


[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn ); %Denne skal byttes ut med "getTestFunctions(prob,X,T,var)"
V(:,1) = A*V(:,1); % denne kodelinjen er teit?
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAN:
% funksjon
% Verifiser inndata
% Initial stuff X
% Få tak i testfunksjoner X--
% en forløkke som løser problemet
% beregne feil X--
% annet

if alg == 1 || alg == 2
    if alg == 1
        alg = @Arnoldi;
    elseif alg == 2
        alg = @SymplecticLanczosMethod;
    end
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
elseif alg == 3
    utdata(1) = 0;
    tic;
    tempVF = 0;
    for i = 1:size(F,1)
        tempVF = tempVF + V(:,i)*F(i,:);
    end
    Utemp = integrate(A,tempVF,k,ht);
    utdata(2) = toc;
    utdata(5) = -1;
end

U = zeros(m^2,k);
Utemp = Utemp + U0*ones(1,k);
U(vec,:) = Utemp(1:(m-2)^2,:);
 if utdata(5) == 0
    tempVF1 = 0;
    for i = 1:size(F,1)
        tempVF1 = tempVF1 + V(:,i)*F(i,:);
    end
    Utemp1 = integrate(A,tempVF1,k,ht);
    U1 = zeros(m^2,k);
    Utemp1 = Utemp1 + U0*ones(1,k);
    U1(vec,:) = Utemp1(1:(m-2)^2,:);
    utdata(5) = max(max(abs(U-U1)));
 end
 
%utdata(3) = getError(U,correctsolution);
utdata(3) = max(max(abs(U-correctsolution)));

if 0
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


%error =
%utdata;
%utdata(3)
%toc;
end