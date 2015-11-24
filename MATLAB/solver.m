
function utdata = solver(m,n,k,eqn,alg,restart,prob,conv,para)
%Rekkefølge inn: m,n,k,eqn,alg,restart,prob,conv,para
%Skriv en programdefinosjon her
%clear
%close all
%%% Initiell data
%tic;
if nargin < 9
    m = 6;
    k = 10;
    n = 4;%2*(m-2)^2;
    restart = 1;
    prob = 1;
    conv = 10^-14;
    para = 4; %%%%% ARG %%%%%%
    eqn = 'maxwell1D';
    alg = 3;
end
%Alt over dette burde være argumenter +
% if ~verifyData(m,n,k,eqn,alg,restart,prob,para)
%     utdata = -ones(1,4);
%     return
% end
%m,n,k,eqn,alg,restart,prob,conv,para
%datastore1 = datastorage(m,n,k,restart,eqn,prob,para,alg,conv,1) % Sjekker om krylov må
%benyttes
%datastore2 = datastorage(m,n,k,restart,eqn,prob,para,alg,conv,2) % Sjekker
%om den er løst av direkte integrasjons metoden


%%% Initsiell data
utdata = zeros(1,5); % Burde legge til forskjellen mellom energi og
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);

[vec,height] = helpvector(m,eqn);

%disk = ht^2/(hs^2);


%%%%%%%%%% TODO %%%%%%%%%%
%%% Skrive bølge, maxwell og varme sammen så mye som mulig %%%
%%% Legge til beskrivelse til ferdige funksjoner %%%
%%% Lage en database av resultater, så beregninger går fortere!
%%% Feilen mellom direkte metode og Krylov metode burde også være med,
%%% alltid!


[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V(:,1) = A*V(:,1);
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
        algo = @Arnoldi;
    elseif alg == 2
        algo = @SymplecticLanczosMethod;
    end
    tic;
    iter = 0;
    Utemp = 0;
    for i = 1:size(F,1)
        [Utemp1,iter1] = KPM(A,V(:,i),F(i,:),n,ht,conv,restart,algo);
        Utemp = Utemp + Utemp1;
        iter = max(iter1,iter);
    end
    utdata(1) = iter;
    utdata(2) = toc;
    U = zeros(height,k);
    Utemp = Utemp + U0*ones(1,k);
    if strcmp(eqn,'maxwell1D')
        U(vec,:) = Utemp(1:length(A)/2-1,:); % OBS: Dette er en dårlig løsning!
    else
        U(vec,:) = Utemp(1:length(A)/2,:);
    end
end


tic;
tempVF = 0;
for i = 1:size(F,1)
    tempVF = tempVF + V(:,i)*F(i,:);
end
Utemp1 = integrate(A,tempVF,k,ht);
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
else
    utdata(1) = 0;
    utdata(2) = Time;
    utdata(5) = -1;
    Utemp = Utemp1;
    U = U1;
end


%utdata(3) = getError(U,correctsolution);
utdata(3) = max(max(abs(U-correctsolution)));
utdata(4) = energy(A,Utemp);
if 1
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