
function utdata = solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,figvar)
% Solves a problem dependant on the indata
%input
% m: number of points in eqch spacial direction X
% n: restart variable (size of orthogonal space) X
% k: number of space in time. X
% eqn: says something about with algorithm to solve
% alg(1,2,3): declares with ortogonalisation method to use X
% integrator(1,2,3): declares with integration method to use X
% restart(0,1): should the method restart or not X
% prob: says smoething about with particular problem to solve
% conv: convergence criterion used in arnoldi og KPM X
% para: currently nothing
%returns:
% utdata:
% utdata(1): number of ioeration performed
% utdata(2): computation time
% utdata(3): error
% utdata(4): energy
% utdata(5): error difference
% utdata(6): energy difference

%%%% ISSUES %%%%%
% Burde alle bildene numeres slik at det er lettere for de forskjellige
% fuksjonene å finne fram det relevante bildet?


%%% Initiell data
if nargin < 10
    m = 20;
    simtime = 10;
    K = 2;
    k = 20;
    n = 6;%2*(m-2)^2;
    restart = 1;
    prob = 1;
    conv = 10^-14;
    para = 4; %%%%% If need be %%%%%%
    eqn = 'wave';
    alg = 2;
    integrator = 3;
    figvar = 1;
end

%%%% lage en funksjon som tar seg av dette? Og alt tilknyttet dette?
% Chose integration method
if integrator == 1
    int = @trapezoidal;
    timestep = 1;
elseif integrator == 2
    int = @forwardeuler;
    timestep = 1;
elseif integrator == 3
    int = @implicitmidpoint;
    k = 2*k;
    timestep = 2;
end

%%% Initsiell data
utdata = zeros(1,6);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,simtime,K*k); ht = T(2)-T(1);
[vec,height,lastrelevant] = helpvector(m,eqn);

% Get problem information
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );

% Se om det går ann å skrive saken nedenfor litt penere!!!!
% Chose solution method and solve
if alg == 1 || alg == 2
    if alg == 1
        algo = @Arnoldi;
    elseif alg == 2
        algo = @SymplecticLanczosMethod;
        n = n/2;
    end
    
    
    tic;
    % Denne algoritmen er noe ROOOT!!!! FIX it!
    iter = 0;
    Utemp = zeros(size(V,1),size(F,2));
    U0temp = U0;
    for j = 1:K
        d = timestep*ceil((K-j)/K); timeinterval = 1+(j-1)*k:j*k+d;
        V(:,1) =  A*U0temp;
        for i = 1:size(F,1)
            [Utemptemp,iter1] = KPM(A,V(:,i),F(i,timeinterval),n,ht,conv,restart,algo,int);
            Utemp(:,timeinterval) = Utemp(:,timeinterval) + Utemptemp;
            iter = max(iter1,iter);
        end
        utdata(1) = max(iter,utdata(1));
        Utemp(:,timeinterval) = U0temp*ones(1,k+d) + Utemp(:,timeinterval);
        U0temp = Utemp(:,j*k+d/max(1,d));
        Utemp(:,j*k+d/max(1,d)) = zeros(size(V,1),1);
    end
    Utemp(:,j*k+d/max(1,d)) = U0temp;
    utdata(2) = toc;
    
    % Det som står under her vil jeg fjerne!!!!!!!!!!%%%%%
    U = zeros(height,K*k/timestep);
    Utemp = Utemp(:,1:timestep:end);
    correctsolution = correctsolution(:,1:timestep:end); T = T(1:timestep:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    U(vec,:) = Utemp(1:lastrelevant,:);
end


% denne kan også gjerne skrives penere
tic;
VF = 0;
for i = 2:size(F,1)
    VF = VF + V(:,i)*F(i,:);
end
Utemp1 = zeros( size(VF) );
U0temp = U0;
for j = 1:K
    d = timestep*ceil((K-j)/K); timeinterval = 1+(j-1)*k:j*k+d;
    tempVF = VF(:,timeinterval) + A*U0temp*ones(1,k+d);
    Utemp1(:,timeinterval) = U0temp*ones(1,k+d) + int(A,tempVF,ht);
    U0temp = Utemp1(:,j*k+d/max(1,d));
end
Time = toc;


% Alt her er nokså teit!!!!!!%%%%%%
U1 = zeros(height,K*k/timestep);
Utemp1 = Utemp1(:,1:timestep:end); k = k/timestep;
if alg == 3
    T = T(1:timestep:end);
    correctsolution = correctsolution(:,1:timestep:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U1(vec,:) = Utemp1(1:lastrelevant,:);

% Plots and outdata
% Det kunne vært en funksjon for plot og utdata?

if alg ~= 3
    utdata(5) = max(max(abs(U-U1)));
    %figure(5);plot(T,max(abs(U-U1)), 'k:.')
    
    %figure(7); plot(T,energy(A,U(vec,:)-U1(vec,:)),'k:.');
    utdata(6) = max(abs(energy(A,U(vec,:)-U1(vec,:))));
else
    utdata(1) = 0;
    utdata(2) = Time;
    utdata(5) = -1;
    utdata(6) = -1;
    U = U1;
end


utdata(3) = max(max(abs(U-correctsolution)));
if figvar
    figure(11); plot(T,max((U-correctsolution)),'k:.');
end

if prob == 1
    utdata(4) = max(abs(energy(A,U(vec,:))));
    if figvar
        figure(2); plot(T,energy(A,U(vec,:)),'k:.');
    end
else
    if figvar
        figure(2); plot(T,energy(A,U(vec,:)-correctsolution(vec,:)),'k:.');
    end
    utdata(4) = abs(max(abs(energy(A,U(vec,:)-correctsolution(vec,:)))));
end

%energy(A,U(vec,:)-correctsolution(vec,:),U0)
%energy(A,U(vec,:),U0)
%energy(A,correctsolution(vec,:),U0)
if figvar
    figure(11);plot(T,max(U-correctsolution),'k:.');
end

% Plot
if 0
    %video(U(m^2+1:end,:),m,0.05,eqn)
    %video(U(1:m^2,:),m,0.05,eqn)
    %video(U1(m^2+1:end,:),m,0.05,eqn)
    %video(U(1:m^2,:)-U1(1:m^2,:),m,0.05,eqn)
    %video(U(m^2+1:end,:)-correctsolution(m^2+1:end,:),m,0.05,eqn)
    %video(U(m^2+1:end,:),m,0.05,eqn)
    %video(correctsolution(m^2+1:end,:)-U(m^2+1:end,:),m,0.05,eqn)
    %video(correctsolution(1:m^2,:)-U(1:m^2,:),m,0.5,eqn)
end
end