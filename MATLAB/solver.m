
function utdata = solver(m,n,k,eqn,alg,integrator,restart,prob,conv,para)
%Rekkefølge inn: m,n,k,eqn,alg,restart,prob,conv,para
%Skriv en programdefinosjon her
%clear
%close all
%%% Initiell data
%tic;
if nargin < 10
    m = 20;
    k = 20;
    n = m;%2*(m-2)^2;
    restart = 0;
    prob = 3;
    conv = 10^-14;
    para = 4; %%%%% If need be %%%%%%
    eqn = 'wave';
    alg = 1;
    integrator = 1;
end

%%% Initsiell data
utdata = zeros(1,6);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);
[vec,height] = helpvector(m,eqn);




%%%%%%%%%% TODO %%%%%%%%%%
%%% Skrive bølge, maxwell og varme sammen så mye som mulig %%%
%%% Legge til beskrivelse til ferdige funksjoner %%%
%%% Lage en database av resultater, så beregninger går fortere!
%%% Feilen mellom direkte metode og Krylov metode burde også være med,
%%% alltid!

% Get problem information
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V(:,1) = A*V(:,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLAN:
% funksjon
% Initial stuff X
% Få tak i testfunksjoner X--
% en forløkke som løser problemet
% beregne feil X--
% annet

if integrator == 1
    int = @trapezoidal;
elseif integrator == 2
    int = @forwardeuler;
elseif integrator == 3
    int = @implicitmidpoint;
end


% Chose solution method and solve
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
    for i = 1:size(F,1)
        [Utemp1,iter1] = KPM(A,V(:,i),F(i,:),n,ht,conv,restart,algo,int);
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
Utemp1 = int(A,tempVF,k,ht);
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
    utdata(6) = abs(energy(A,Utemp-Utemp1,T));
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
utdata(4) = energy(A,Utemp,T);
figure(5);plot(T,max(U-U1), 'k:.')
% Plot
if 1
    %V = zeros(m^2,k);
    %V(vec,:) = Utemp((m-2)^2+1:end,:);
    %V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);
    %video(U,m,k,0.05,eqn)
    video(U1-U,m,k,0.05,eqn)
    %video(V,m,k,0.05)
    %video(correctsolution,m,k,0.05,eqn)
    %video(U-correctsolution,m,k,0.05,eqn)
    %energy(Jtilde*Atilde,Utemp);
end
end