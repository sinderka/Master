function utdata = energyTest(m,n,simtime,k,eqn,integrator,restart,prob,conv,figvar,solveexpm,intesolve,SLMint, save)
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
% utdata(3): error
% utdata(4): energy
% utdata(5): error difference
% utdata(6): energy difference
% Initsiell
%close all
if nargin < 9
    m = 30;
    simtime = 10;
    %K = 2;
    k = 500;
    n = 4;%2*(m-2)^2;
    restart = 0 ;
    prob = 2;
    conv = 10^-6;
    para = 4; %%%%% If need be %%%%%%
    eqn = 'random';
    integrator = 1;
    figvar = 1;
    solveexpm = 0;
    intesolve = 1;
    PMint = 3;
end

if integrator == 1
    int = @trapezoidal;
elseif integrator == 2
    int = @forwardeuler;
end

utdata = zeros(1,6);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,simtime,k); ht = T(2)-T(1);
%[vec,height,lastrelevant] = helpvector(m,eqn);

% Get problem information
[A] = getMatrix( m , hs, eqn );
[U0,V,F,correctsolution] = getTestFunctions( prob,X,T,eqn );
V = A*U0;

% Chose solution method and solve
%U = zeros(2*m^2,k);

hamiltonian(A);

%energy1 = 0; energy2 = 0; iter = 0;
tic;
%for i = 1:size(F,1)
[Ukpm,~,~,~,~,~,~,~] = PM(A,V,ones(1,k),n,ht,conv,restart,int,0,PMint,@Arnoldi);
[Uslm,~,~,~,~,~,~,~] = PM(A,V,ones(1,k),2*n,ht,conv,restart,int,0,PMint,@SLM);
 %   [Utemp,iter1,energy1t,energy2t] = PM(A,V,T,n/2,conv,restart,ht,figvar,int,SLMint,correctsolution(vec,:)-U0*ones(1,k));
 %   [Utemp,iter1,energy1t,energy2t] = PM(A,V,T,n/2,conv,restart,ht,figvar,int,SLMint,correctsolution(vec,:)-U0*ones(1,k));
 %   U(vec,:) = U(vec,:) + Utemp;
 %   energy1 = energy1 + energy1t; energy2 = energy2 + energy2t; iter = max(iter1,iter);
%end

Ukpm = Ukpm + U0*ones(1,k);
Uslm = Uslm + U0*ones(1,k);
utdata(2) = toc;

%utdata(9) = energy1; utdata(10) = energy2;
%utdata(1) = iter;

    %expmsolution = zeros(2*m^2,k);
    expmsolution = expintegrate(A,U0,T) + U0*ones(1,k);
    %utdata(5) = max(getError(U,expmsolution));
    %utdata(6) = max(abs(getEnergy(A,U(vec,:))-getEnergy(A,expmsolution(vec,:))));


if intesolve && 0
    intesolution = zeros(2*m^2,k);
    intesolution(vec,:) = int(A,V(:,1)*ones(1,k),ht);% + U0*ones(1,k);
    utdata(7) = max(getError(U,intesolution));
    utdata(8) = max(abs(getEnergy(A,U(vec,:))-getEnergy(A,intesolution(vec,:))));
end
for i = 1:k
    errslm(i) = norm(expmsolution(:,i)-Uslm(:,i))/norm(expmsolution(:,i));
    errkpm(i) = norm(expmsolution(:,i)-Ukpm(:,i))/norm(expmsolution(:,i));
end
figure(1)
plot(errslm)
hold on
plot(errkpm)
legend('SLM','KPM')

hold off
for i = 1:k
    
    eneslm(i) = 1/2*(Uslm(:,i)'*A*Uslm(:,i));
    enekpm(i) = 1/2*(Ukpm(:,i)'*A*Ukpm(:,i));
    
    %eneslm(i) = 1/2*(expmsolution(:,i)'*A*expmsolution(:,i)-Uslm(:,i)'*A*Uslm(:,i));
    %enekpm(i) = 1/2*(expmsolution(:,i)'*A*expmsolution(:,i)-Ukpm(:,i)'*A*Ukpm(:,i));
end
figure(2)
plot(eneslm)
hold on
plot(enekpm)
legend('SLM','KPM')



if figvar&&0
    figure(111); plot(T,getError(U(vec,:),correctsolution(vec,:)),'k:.')
    %figure(13); plot(T,getEnergy(A,U(vec,:)),'k:.');
    figure(12); plot(T,getEnergy(A,U(vec,:),V(:,1)),'k:.'); %hold on; plot(T,T*1e-12,'k-'); hold off;
    if solveexpm
        figure(15); plot(T,getError(U,expmsolution),'k:.')
        figure(17); loglog(T,getEnergy(A,expmsolution(vec,:)),'k:.'); %hold on; plot(T,T.^2*1e-23,'k-'); hold off;
    end
    if intesolve
        figure(18); plot(T,getError(U,intesolution),'k:.')
        figure(19); plot(T,getEnergy(A,intesolution(vec,:),V(:,1)),'k:.'); %hold on; plot(T,T.^2*1e-23,'k-'); hold off;
        figure(20); plot(T,getError(intesolution(vec,:),correctsolution(vec,:)),'k:.')
    end
    
end
utdata(3) = max(getError(U,correctsolution));
utdata(4) = max(getEnergy(A,U(vec,:)));
a = 2;
%figure(31);loglog(T,getEnergy(A,U(vec,:)),'k:+'); hold on; %plot(T,T.^0*1e-13,'k-'); hold off;
% %getLabels(1,m,n,simtime,1,k,'wave',2,integrator,restart,1,conv,para,4) ; saveit(strcat('energytest1',num2str(n),num2str(SLMint),num2str(simtime)),'T_s','en_1');
%
% figure(32); loglog(T,abs(energy(A,expmsolution(vec,:))),'ko:');% hold on; plot(T,T.^0*1e-12,'k-'); hold off;
% %getLabels(1,m,n,simtime,1,k,'wave',2,integrator,restart,1,conv,para,4) ; saveit(strcat('energytest2',num2str(n),num2str(SLMint),num2str(simtime)),'T_s','en_1');
%
% figure(33);
% loglog(T,abs(energy(A,intesolution(vec,:))),'kx:');% hold on; plot(T,T.^0*1e-13,'k-'); hold off;
% %[~, ~,~,additionalInfo] = getLabels(1,m,n,simtime,1,k,'wave',2,integrator,restart,1,conv,1,4); legend('SLM','EXPm','intmeth'); title(additionalInfo);
% if save
%     saveit(strcat('energytest',num2str(n),num2str(SLMint),num2str(simtime),num2str(integrator)),'T_s','en_1');
% end
% hold off;

if 0
    %video(U(1:m^2,:),m,0.05,eqn)
    %pause
    video(U(1:m^2,:)-correctsolution(1:m^2,:),m,0.5,eqn)
    
    %video(expmsolution(1:m^2,:),m,0.05,eqn)
    %video(expmsolution(1:m^2,:)-correctsolution(1:m^2,:),m,0.05,eqn)
    
    %video(intesolution(1:m^2,:),m,0.05,eqn)
    video(intesolution(1:m^2,:)-correctsolution(1:m^2,:),m,0.5,eqn)
    
    %video(correctsolution(1:m^2,:),m,0.05,eqn)
    
    %video(U(1:m^2,:)-expmsolution(1:m^2,:),m,0.05,eqn)
end
end

function [U,iter,energy1,energy2] = KPMloc(A,v,T,n,conv,restart,ht,figvar,int,SLMint,correctsolution)
%Indata
% A: mxm matrix
% v: m vector
% F: k row of timedependant function
% n: real number 0<n<=m
% ht: stepsize in time
% conv: convergence criterion
% restart: A boolean value
% alg: an ortogonalisation algorithm (Arnoldi or SLM)
% int: an integration method (trapezoidal rule)
%outdata
% U: Solution to problem du/dt = Au+v*F
% iter: number of restarts preformed
l = size(A,1);
k = length(T);
if max(abs(v)) == 0
    U = sparse(l,k);
    iter = 0;
    energy1 = 0;
    energy2 = 0;
    return
end
U = zeros(l,k);
iter = 1;
%v0 = v;
h = norm(v,2);

[Vn,Hn,vnext,hnext] = SLM(A,v,n,conv);
hamiltonian(Hn); symplectic(Vn); eigenschaft(A,Vn,Hn,vnext,hnext);

if SLMint == 1
    Zn = int(Hn,[h*ones(1,k);sparse(length(Hn)-1,k)],ht);
elseif SLMint == 2
    Zn = expintegrate(Hn,Hn\[h;sparse(length(Hn)-1,1)],0:ht:ht*(k-1));
elseif SLMint == 3
    Zn = real(myexpm(full(Hn),Hn\[h;sparse(length(Hn)-1,1)],0:ht:ht*(k-1)));
end


figure(231);loglog(0:ht:ht*(k-1),getEnergy(Hn,Zn,[h;sparse(length(Hn)-1,1)]),'k:.')
%figure(234);loglog(0:ht:ht*(k-1),getEnergy(Hn,Zn),'k:.')

%ns = int(A,v*ones(1,k),ht); figure(232);loglog(0:ht:ht*(k-1),getEnergy(A,ns,v),'k:.')

ns = Vn*Zn;

assumtion(A,v,Zn,correctsolution)

Vn0 = Vn; Zn0 = Zn; Hn0 = Hn; v0 = v;
U = U + ns;

figure(232);loglog(0:ht:ht*(k-1),getEnergy(A,U,v),'k:.')
% epsilon = correctsolution-U;
% something = zeros(1,k);
% btilde = [norm(v0,2)*ones(1,k);sparse(length(Hn)-1,k)];
% J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
% for kk = 1:k
%     something(:,kk) = epsilon(:,kk)'*J*Vn*(Hn*Zn(:,kk) +btilde(:,kk) );
% end
% figure(312); plot(something)

            %energy1 = max(abs(energyBIG(A,Zn,v,h,ht,figvar,int)));
            %energy2 = max(abs(energySMALL(Hn,Vn,Zn,v,h,ht,figvar,int)));

diff = hnext;
if restart
    
    %invJ = [sparse(n,n),-speye(n);speye(n),sparse(n,n)];
    %J = [sparse(l/2,l/2),speye(l/2);-speye(l/2),sparse(l/2,l/2)];
    
    while diff > conv
        h = hnext; v = vnext;
        [Vn,Hn,vnext,hnext] = SLM(A,v,n,conv);
        hamiltonian(Hn); symplectic(Vn); eigenschaft(A,Vn,Hn,vnext,hnext);
        
        Zn = int(Hn,[h*ones(1,k);sparse(length(Hn)-1,k)],ht);
        %Zn = int(Hn,F,ht);
        ns =  Vn*Zn;
        diff = max(max(abs(ns)));
        U = U + ns;
        iter = iter+1;
        if iter == 2
            luliproof(A,v0,h*v,Vn0,Vn,Hn0,Hn,Zn0,Zn);
            energy1 = max(abs(energyBIG(A,Zn,v,h,ht,figvar,int)));
            energy2 = max(abs(energySMALL(Hn,Vn,Zn,v,h,ht,figvar,int)));
            stuff = energySMALL(Hn,Vn,Zn,v,h,ht,figvar,int); things = energyBIG(A,Zn,v,h,ht,figvar,int);
            figure(21983);plot(stuff- things,'k:.');
            %energy2 = max(abs(energySMALL(Hn,Vn,Zn,v,h,ht,figvar,int)));
        end
    end

end
if ~exist('energy1','var')
    energy1 = -1;
    energy2 = -1;
end

%stuff = energySMALL(Hn,Vn,Zn,v,h,ht,figvar,int); things = energyBIG(A,Zn,v,h,ht,figvar,int);
%figure(21983);plot(stuff- things);

% F = invJ*Vn'*J*hnext*vnext*Zn(end,:);
% 
% delta = int(Hn,F,ht);
% blah =  Vn0*Zn0 + Vn*delta ;
% if figvar
%     figure(16);plot(getEnergy(A,blah,v0),'k:.')
% end

end

function hamiltonian(A)
m = length(A)/2;
J = [sparse(m,m),speye(m); -speye(m),sparse(m,m)];
if isequal((J*A)',J*A);
    display('A er Hamiltonsk');
else
    display('A er ikke Hamiltonsk');
end
end

function symplectic(A)
[m,n] = size(A);
Jm = [sparse(m/2,m/2),speye(m/2); -speye(m/2),sparse(m/2,m/2)];
Jn = [sparse(n/2,n/2),speye(n/2); -speye(n/2),sparse(n/2,n/2)];
if isequal((A'*Jm*A)',Jn)%max(max(abs((A'*Jm*A)'+Jn))) < 1e-15
    display('S_n er symplektisk');
else%if max(max(abs((A'*Jm*A)'+Jn))) < 1e-15
    num = max(max(abs((A'*Jm*A)'+Jn)));
    fprintf('Den maksimale forskjellen mellom (S^T*Jm*S)^T og Jn er %d\n',num);
end
end

function eigenschaft(A,Sn,Hn,vnext,hnext)
e_n = zeros(1,length(Hn)); e_n(end) = 1;
%e_1 = zeros(1,length(Hn)); e_n(end) = 1;

if A*Sn-Sn*Hn-hnext*vnext*e_n < 10^-10
    display('SLM virker')
else
    display('SLM virker ikke')
    
end
end

function U = locexpm(A,b,T)
U = zeros(length(A),length(T));
[V,D] = eig(A);
for i = 1:length(T)
    U(:,i) = V*diag(exp(diag(D*T(i))))/V * b - b;
end
end

function saveit(name,xlab,ylab)
xlabel(xlab)
ylabel(ylab)
h = set(findall(gcf,'-property','FontSize'), 'Fontsize',12);
set(h,'Location','Best');
pause(0.5)
drawnow
pause(0.5)
location = strcat('/home/shomeb/s/sindreka/Master/MATLAB/fig/',char(name));
saveas(gcf,location,'fig');
saveas(gcf,location,'jpeg');

end


function luliproof(A,b,b2,Sn1,Sn2,Hn1,Hn2,z,delta)
%z, Sn1 er for å løse forste gangen
%Sn2 og delta er fra å løse senere ganger(andre gangen)

%J og A er greie!

l = size(A,1); [n,k] = size(z);
J = [sparse(l/2,l/2),speye(l/2,l/2);-speye(l/2,l/2),sparse(l/2,l/2)];
invJ = [sparse(n/2,n/2),-speye(n/2,n/2);speye(n/2,n/2),sparse(n/2,n/2)];
Jn = [sparse(n/2,n/2),speye(n/2,n/2);-speye(n/2,n/2),sparse(n/2,n/2)];
btilde = zeros(n,1);btilde(1) = norm(b,2);%invJ*Sn1'*J*b;
zdot = Hn1*z+btilde*ones(1,k);

e2n = zeros(n,1); e2n(end) = 1;
%e2n = zeros(n,1); e2n(1) = 1;

Sn1ting = (Sn1*z+Sn2*delta);
var1 = zeros(1,k); var2 = zeros(1,k); var3 = zeros(1,k); var4 = zeros(1,k); var6 = zeros(1,k);
for kk = 1:k
    var1(kk) = 1/2*Sn1ting(:,kk)'*J*A*Sn1ting(:,kk) + Sn1ting(:,kk)'*J*b;%*ones(1,k);
    var2(kk) = (Sn2*delta(:,kk))'*J*Sn1*zdot(:,kk);
    var6(kk) = 1/2*delta(:,kk)'*Jn*Hn2*delta(:,kk) + delta(:,kk)'*Sn2'*J*b2*e2n'*z(:,kk);
    var3(kk) = 1/2*delta(:,kk)'*Jn*Hn2*delta(:,kk) + (Sn2*delta(:,kk))'*J*(Sn1*Hn1+b2*e2n')*z(:,kk) + delta(:,kk)'*Sn2'*J*b;
    var4(kk) = 1/2*delta(:,kk)'*Jn*Hn2*delta(:,kk) + delta(:,kk)'*Sn2'*J*b2*e2n'*z(:,kk);
end
if max(max(abs(var1-var2))) < 1e-10
    display('!!!!!!!!!!!!!!LuLi ting fungerer!!!!!!!!!!!!!!!!')
else
    display('!!!!!!!!!!!!!!!!!!!!!!!!!Sjekk for feil!!!!!!!!!')
end


end

function assumtion(A,v,Zn,correctsolution)

[Vn,Hn,vnext,hnext] = SLM(A,v,size(A,2)/2,1);

Z = Vn'*correctsolution;

for i = 1:size(Zn,1)
    if max(abs(Zn(i,:) - Z(i,:))) < 1e-10
        fprintf('Det gikk bra for i = %d \n',i)
    else
        fprintf('Forskjellen for i = %d er: %d\n',i,max(abs(Zn(i,:) - Z(i,:))))
    end
end




end