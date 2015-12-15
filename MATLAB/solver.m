
function utdata = solver(m,n,k,eqn,alg,integrator,restart,prob,conv,para)
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

%%% Initiell data
if nargin < 10
    m = 20;
    k = 20;
    n = 6;%2*(m-2)^2;
    restart = 0;
    prob = 1;
    conv = 10^-14;
    para = 4; %%%%% If need be %%%%%%
    eqn = 'wave';
    alg = 3;
    integrator = 1;
end

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
T = linspace(0,1,k);ht = T(2)-T(1);
[vec,height,lastrelevant] = helpvector(m,eqn);

% Get problem information
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V(:,1) = A*V(:,1);




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
        [Utemptemp,iter1] = KPM(A,V(:,i),F(i,:),n,ht,conv,restart,algo,int);
        Utemp = Utemp + Utemptemp;
        iter = max(iter1,iter);
    end
    utdata(1) = iter;
    utdata(2) = toc;
    
    
    % Det som står under her vil jeg fjerne!!!!!!!!!!%%%%%
    U = zeros(height,k/timestep);
    Utemp = Utemp(:,1:timestep:end);
    correctsolution = correctsolution(:,1:timestep:end); T = T(1:timestep:end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Utemp = Utemp + U0*ones(1,k/timestep); % Pluss noe her
    U(vec,:) = Utemp(1:lastrelevant,:);
end


tic;
tempVF = 0;
for i = 1:size(F,1)
    tempVF = tempVF + V(:,i)*F(i,:);
end
Utemp1 = int(A,tempVF,ht);
Time = toc;


Utemp1 = Utemp1 + U0*ones(1,k);

% Alt her er nokså teit!!!!!!%%%%%%
U1 = zeros(height,k/timestep);
Utemp1 = Utemp1(:,1:timestep:end); k = k/timestep;
if alg == 3
    T = T(1:timestep:end);
    correctsolution = correctsolution(:,1:timestep:end);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U1(vec,:) = Utemp1(1:lastrelevant,:);

% Plots and outdata

if alg ~= 3
    utdata(5) = max(max(abs(U-U1)));
    figure(5);plot(T,max(abs(U-U1)), 'k:.')
    
    figure(7); plot(T,energy(A,Utemp-Utemp1,U0),'k:.');
    utdata(6) = max(abs(energy(A,Utemp-Utemp1,U0)));
else
    utdata(1) = 0;
    utdata(2) = Time;
    utdata(5) = -1;
    utdata(6) = -1;
    Utemp = Utemp1;
    U = U1;
end


utdata(3) = max(max(abs(U-correctsolution)));
figure(11); plot(T,max(abs(U-correctsolution)),'k:.');

figure(2); plot(T,energy(A,Utemp,U0),'k:.');
utdata(4) = max(abs(energy(A,Utemp,U0)));


figure(11);plot(T,max(abs(U-correctsolution)),'k:.');

% Plot
if 0
    %V = zeros(m^2,k);
    %V(vec,:) = Utemp((m-2)^2+1:end,:);
    %V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);
    video(U,m,k,0.05,eqn)
    %video(U1-U,m,k,0.05,eqn)
    %video(V,m,k,0.05)
    %video(correctsolution,m,k,0.05,eqn)
    %video(U-correctsolution,m,k,0.05,eqn)
    %energy(Jtilde*Atilde,Utemp);
end
end