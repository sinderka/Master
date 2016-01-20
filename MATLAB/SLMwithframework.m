function utdata = SLMwithframework(m,n,simtime,k,integrator,restart,conv,figvar,solveexpm)

% The relevnat functions are:
% this function
% intloc
% PMloc
% trapezoidalLoc
% forwardeulerLoc

% This function shopuld work now, you can add/change problem in get
% Testproblems and getMatrix.

% Solves a problem dependant on the indata
%input
% m: number of points in eqch spacial direction X
% n: restart variable (size of orthogonal space) X
% k: number of space in time. X
% eqn: says something about with algorithm to solve
% alg(1,2,3): declares with ortogonalisation method to use X
% integrator(1,2,3): declares with integration method to use X
% restart(0,1): should the method restart or not X
% prob: says something about with particular problem to solve
% conv: convergence criterion used in arnoldi og KPM X
% para: currently nothing
%returns:
% utdata:
% utdata(1): number of iteration performed
% utdata(2): computation time
% utdata(3): global error
% utdata(4): global energy
% utdata(5): error compared with expmsolution
% utdata(6): energy compared with expmsolution
% utdata(7): error compared with just integration
% utdata(8): energy compared with just integration
% Initsiell
if nargin < 9
    m = 20; % Number of points in each spacial direction
    simtime = 1; % number of seconds simulated
    k = 20; % Number of steps in time
    n = 2; % restart variable 2,4,6,...,2(m-2)^2
    restart = 1; % Enable restart ( must be either 0 or 1)
    conv = 10^-10; % convergence criterion or tolerance
    integrator = 1; % Integration method (1 == trapezoidal, 2 == forwardeuler)
    figvar = 1; % Make figures
    solveexpm = 1; % solve problem with exponential function (may be time consuming)
    solveinte = 1; % solve problem with integration function
end

% Choose integration method
if integrator == 1
    int = @trapezoidalLoc; % Loc means it is a local function
elseif integrator == 2
    int = @forwardeulerLoc;
end

% initial data
utdata = zeros(1,8); utdata(5:end) = -1*ones(1,4);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,simtime,k); ht = T(2)-T(1);
[vec] = helpvectorLoc(m);
[A] = getMatrixLoc( m , hs );
[U0,V,correctsolution] = getTestFunctionsLoc(X,T);
V = A*V;

% Use projection method to solve the problem
U = zeros(2*m^2,k);
iter = 0;
tic;
[Utemp,iter1] = PMloc(A,V,T,n/2,conv,restart,ht,int);
U(vec,:) = U(vec,:) + Utemp;
iter = max(iter1,iter);
U(vec,:) = U(vec,:) + U0*ones(1,k);

% save outdata
utdata(1) = iter;
utdata(2) = toc;
utdata(3) = max(getError(U,correctsolution));
utdata(4) = max(abs(energy(A,U(vec,:))));

% Solve with exponential function
if solveexpm
    expmsolution = zeros(2*m^2,k);
    expmsolution(vec,:) = intloc (A,U0,T) + U0*ones(1,k);
    utdata(5) = max(getError(U,expmsolution));
    utdata(6) = max(abs(energy(A,U(vec,:)-expmsolution(vec,:))));
end

% Solve with integration function
if solveinte
    intesolution = zeros(2*m^2,k);
    intesolution(vec,:) = int(A,V*ones(1,k),ht) + U0*ones(1,k);
    utdata(7) = max(getError(U,intesolution));
    utdata(8) = max(abs(energy(A,U(vec,:)-intesolution(vec,:))));
end


% Plot figures
if figvar
    figure(2); plot(T,getError(U,correctsolution),'k:.')
    figure(3); plot(T,abs(energy(A,U(vec,:))),'k:.');
    if solveexpm
        figure(4); plot(T,getError(U,expmsolution),'k:.')
        figure(5); plot(T,energy(A,U(vec,:)-expmsolution(vec,:)),'k:.');
    end
    if solveinte
        figure(6); plot(T,getError(U,intesolution),'k:.')
        figure(7); plot(T,energy(A,U(vec,:)-intesolution(vec,:)),'k:.');
    end
    
end

if 1 % See the video of the solution
    videoLoc(U(1:m^2,:),m,0.05)
    videoLoc(U(1:m^2,:)-correctsolution(1:m^2,:),m,0.05)
    
    videoLoc(expmsolution(1:m^2,:),m,0.05)
    videoLoc(expmsolution(1:m^2,:)-correctsolution(1:m^2,:),m,0.05)
    
    videoLoc(intesolution(1:m^2,:),m,0.05)
    videoLoc(intesolution(1:m^2,:)-correctsolution(1:m^2,:),m,0.05)
    
    
    videoLoc(correctsolution(1:m^2,:),m,0.05)
    
end
end

function [U,iter] = PMloc(A,v,T,n,conv,restart,ht,int)
%Indata
% A: mxm matrix
% v: m vector
% T: A row of the relenvant time points
% n: real number 0<n<=m
% ht: stepsize in time
% conv: convergence criterion
% restart: A boolean value
% int: an integration method (trapezoidal rule or forward euler)
%outdata
% U: Solution to problem du/dt = Au+v*F
% iter: number of restarts preformed
l = size(A,1);
k = length(T);
if max(abs(v)) == 0
    U = sparse(l,k);
    iter = 0;
    return
end
U = zeros(l,k);
iter = 1;

[Vn,Hn,vnext,hnext] = SymplecticLanczosMethodLoc(A,v,n,conv);

% Initialising som matrices and vectors
invJ = [sparse(n,n),-speye(n);speye(n),sparse(n,n)];
J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
F = invJ*Vn'*J*v;

% Choose a method of integration
%Zn = trapezoidal(Hn,F*ones(1,k),ht);
Zn = intloc(Hn,Hn\F,T);

% ns = next step.
ns = Vn*Zn;
U = U + ns;
diff = hnext;


if restart
    
    % Initialising som matrices
    invJ = [sparse(n,n),-speye(n);speye(n),sparse(n,n)];
    J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
    
    while diff > conv
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = SymplecticLanczosMethodLoc(A,v,n,conv);
        F = invJ*Vn'*J*h*v*Zn(end,:);
        Zn = int(Hn,F,ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
    end
end
end

function U = intloc(A,U0,T)
%solves the problem du/dt = Au+U0, u(0) = 0
%indata
% A: mxm matrix
% u0: initial condition
% k: number of points in time
% ht: step size in time
%outdata
% U: the solution

U = zeros(length(A),length(T));
for i = 1:length(T)
    U(:,i) = expm(A*T(i))*U0-U0;
end
end

function U = trapezoidalLoc( A,F,ht )
%solves the problem du/dt = Au+F
%indata
% A: mxm matrix
% F: mxk matrix
% ht: step size in time
%outdata
% U: the solution
k = size(F,2);
n = size(A,1);
U = zeros(n,k);
mat = speye(n)-A*ht/2;

for i = 2:k
    U(:,i) = mat\(U(:,i-1) + ht/2*A*U(:,i-1)+ht/2*(F(:,i)+F(:,i-1)));
end
end

function U = forwardeulerLoc( A,F,ht )
%solves the problem du/dt = A*u+F
%indata
% A: mxm matrix
% F: mxk matrix
% ht: step size in time
%outdata
% U: the solution
k = size(F,2);
n = size(A,1);
U = zeros(n,k);
for i = 2:k
    U(:,i) = U(:,i-1) + ht*( A*U(:,i-1) + F(:,i-1));
end
end

function [S,Htilde,Vend,xiend] = SymplecticLanczosMethodLoc(H,v,var,~)
% Ortohonalisation method
%Input
% H: a hamiltonian matrix mxm matrix
% v: a m vector
% var: half the size of the resulting orthogonal system
%Returns
% S: a 2*var x m system of orthogonal vectors
% Htilde: a 2*var x 2*var matrix
% Vend: residual vector
% xiend: size of residual vector

n = length(H)/2;
J = [sparse(n,n),speye(n);-speye(n),sparse(n,n)];
delta = zeros(var,1);
beta = zeros(var,1);
xi = zeros(var+1,1);
nu = zeros(var,1);

V = zeros(2*n,var+1);
W = zeros(2*n,var);


xi(2) = norm(v,2);

V(:,2) = 1/xi(2)*v;

for m = 1:1:var
    % Computing v
    v = H*V(:,m+1);
    % Computing delta
    delta(m) = V(:,m+1)'*v;
    % Computing Wm
    wtilde = v-delta(m)*V(:,m+1);
    
    nu(m) = V(:,m+1)'*J*v;
    
    W(:,m) = 1/nu(m)*wtilde;
    W(:,m) = W(:,m)+[V(:,2:m),W(:,1:m-1)]*[sparse(m-1,m-1),speye(m-1);-speye(m-1),sparse(m-1,m-1)]*[V(:,2:m),W(:,1:m-1)]'*J*W(:,m);
    % Computing w
    w = H*W(:,m);
    
    % Computing beta
    beta(m) = -W(:,m)'*J*w;
    %Computing Wm+1
    vmtilde = w-xi(m+1)*V(:,m)-beta(m)*V(:,m+1)+delta(m)*W(:,m);
    xi(m+2) = norm(vmtilde,2);
    V(:,m+2) = 1/xi(m+2)*vmtilde;
    V(:,m+2) = V(:,m+2)+[V(:,2:m+1),W(:,1:m)]*[sparse(m,m),speye(m);-speye(m),sparse(m,m)]*[V(:,2:m+1),W(:,1:m)]'*J*V(:,m+2);
end

S = [V(:,2:end-1),W];

% If var==1 tridiag does not work, therefore it needs to be split in two
% cases
if var > 1
    Htilde = [sparse(1:var,1:var,delta,var,var),gallery('tridiag',xi(3:end-1),beta,xi(3:end-1));
        sparse(1:var,1:var,nu,var,var), sparse(1:var,1:var,-delta,var,var)];
else
    Htilde = [delta,beta;
        nu, -delta];
end


Vend = V(:,end); xiend = xi(end);
end

function [Ustart, V,correctsolution] = getTestFunctionsLoc(X,T)
%Takes 4 agruments;
% X: a set of points in spacial direction
% T: a set of points in time
%Returns 4 matrices
% Ustart: The initial value of the testproblem
% F: A list of vectors depending on time
% V: A list of vectors to generate the Krylov space
% correctsolution: The correct solution.

m = length(X); k = length(T);
[vec] = helpvectorLoc(m);

sol1 = @(t,x,y)  sin(1*pi*x)*sin(1*pi*y)*cos(sqrt(2)*pi*t);
sol2 = @(t,x,y) -sin(1*pi*x)*sin(1*pi*y)*sqrt(2)*pi*sin(sqrt(2)*pi*t);
u0 = @(x,y) sin(1*pi*x)*sin(1*pi*y);
v0 = @(x,y) 0;

Ustart = [getInitial(u0,X);getInitial(v0,X)];
V =  Ustart;
correctsolution = getSolution(X,T,sol1,sol2);

    function correctsolution = getSolution(X,T,sol1,sol2)
        
        correctsolution = zeros(2*m^2,k);
        for j = 1:k
            for i = 1:m
                for l = 1:m
                    correctsolution(l+(i-1)*m,j) = sol1(T(j),X(i),X(l));
                    correctsolution(m^2+l+(i-1)*m,j)=sol2(T(j),X(i),X(l));
                end
            end
        end
        
    end
    function U0 =  getInitial(u0,X)
        U0 = zeros(m^2,1);
        for i = 1:m
            for j = 1:m
                U0(i+(j-1)*m,1) = u0(X(j),X(i));
            end
        end
        U0 = U0(vec(1:length(vec)/2));
    end
end

function [A] = getMatrixLoc( m , hs )
% returns a matrix corresponding to the discretised wave equation
%input:
% m: number of points in space
% hs: stempelngth in space
%returns:
% A 2*(m-2)^2 matrix

A = 1/hs^2*gallery('poisson', m-2);
A = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-A,sparse((m-2)^2,(m-2)^2)];

end

function [v] = helpvectorLoc(m)
% Returns a list of indecis corresponding to non-edge points in a 2 D
% system.
%input:
% m: number of points in eqch spacial direction
%Returns:
% v: a vector og size 2*(m-2)^2 corresponding to non edge points.


v = zeros((m-2)^2,1);
for qq = 0:(m-2)-1
    v(qq*(m-2)+1:qq*(m-2)+m-2) = (qq+1)*m+2:m-1 +(1+ qq)*m;
end
v = [v;v+m^2];
end

function videoLoc(U,m,T)
figure(1)
if max(max(U)) == min(min(U))
    display('Something went wrong with the video!')
    return
end
k = size(U,2);

for i = 1:k
    mesh(reshape(U(:,i),m,m))
    caxis([min(min(U)),max(max(U))])
    axis([0,m,0,m,min(min(U)),max(max(U))])
    drawnow
    pause(T)
end

end







