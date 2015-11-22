function utdata = solver(m,n,k,eqn,alg,restart,prob,conv,para)
%Rekkefølge inn: m,n,k,eqn,alg,restart,prob,conv,para
%Skriv en programdefinosjon her
%clear
%close all
%%% Initiell data
if nargin < 9
    m = 6;
    k = 8;
    n = 6;%2*(m-2)^2;
    restart = 1;
    prob = 1;
    conv = 10^-5;
    para = 4; %%%%% ARG %%%%%%
    eqn = 'wave';
    alg = 3;
end

% if ~verifyData(m,n,k,eqn,alg,restart,prob,para)
%     utdata = -ones(1,4);
%     return
% end

%%% Initsiell data
%utdata = zeros(1,5); % Burde legge til forskjellen mellom energi og
X = linspace(0,1,m); hs = X(2)-X(1);
T = linspace(0,1,k); ht = T(2)-T(1);

vec = helpvector(m,eqn);

%disk = ht^2/(hs^2);


%%%%%%%%%% TODO %%%%%%%%%%
%%% Skrive bølge, maxwell og varme sammen så mye som mulig %%%
%%% Legge til beskrivelse til ferdige funksjoner %%%
%%% Lage en database av resultater, så beregninger går fortere!
%%% Gjøre todo listene på de andre programmene
%%%

% Maxwell skalerer feil med tid OG rom!
% Get crackin'!


[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn,hs );
%V(:,1) = A*V(:,1); V(:,2) = A*V(:,2);

% [bool1,utdata] = datastorage(m,n,k,restart,eqn,prob,para,alg,conv);
% [bool2,U1] = datastorage(m,n,k,restart,eqn,prob,para,alg,conv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAN:
% funksjon
% Verifiser inndata
% Initial stuff X
% Få tak i testfunksjoner X--
% en forløkke som løser problemet
% beregne feil X--
% annet
%[bool1,utdata] = datastorage(m,n,k,restart,eqn,prob,para,alg,conv);
% Enten bruke denne, eller laste den fra et annet sted
if (alg == 1 || alg == 2)% && ~bool1
    if alg == 1
        algtemp = @Arnoldi;
    elseif alg == 2
        algtemp = @SymplecticLanczosMethod;
    end
    tic;
    iter = 0;
    Utemp = 0;
    for i = 1:size(F,1)
        [Utemp1,iter1] = KPM(A,V(:,i),F(i,:),n,ht,conv,restart,algtemp);
        Utemp = Utemp + Utemp1;
        iter = max(iter1,iter);
    end
    utdata(1) = iter;
    utdata(2) = toc;
end

% enten bruke denne, eller laste den fra et annet sted
iter = 0;
tic;
tempVF = 0;
for i = 1:size(F,1)
    tempVF = tempVF + V(:,i)*F(i,:);
end
Utemp1 = integrate(A,tempVF,k,ht);
time = toc;
if alg == 3
    Utemp = Utemp1;
    utdata(1) = iter;
    utdata(2) = time;
    utdata(5) = -1;
    
else
    utdata(5) = max(max(abs(Utemp-Utemp1)));
end

if strcmp(eqn,'maxwell1D')
    he = m;
else
    he = m^2;
end


%U = zeros(he,k);
%Utemp = Utemp + U0*ones(1,k);
%U(vec,:) = Utemp(1:m-2,:);

U = zeros(he,k);
Utemp = Utemp + U0*ones(1,k);
U(vec,:) = Utemp(1:length(A)/2,:);

%utdata(3) = getError(U,correctsolution);
utdata(3) = max(max(abs(U-correctsolution)));
utdata(4) = energy(A,Utemp);


if 1
    %V = zeros(m,k);
    %V(:,:) = Utemp(m-1:end,:);
    %V(1,:) = F(2,:);
    %V(end,:) = -F(2,:);
    %V(:,:) = V(:,:) + U0(m-1:end)*ones(1,k);
    video(U,m,k,0.05,eqn)
    %video(V,m,k,0.05,eqn)
    %video(correctsolution,m,k,0.05,eqn)
    %video(U-correctsolution,m,k,0.05)
    %energy(Jtilde*Atilde,Utemp);
end

% if ~bool1
%     savedatastorage(  m,n,k,restart,eqn,prob,para,alg,conv ,utdata )
% end
end