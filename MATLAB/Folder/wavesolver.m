%Skriv en programdefinosjon her
clear
close all
%%% Initsiell data
tic;
m = 4000;
k = 10;
n = 4;%2*(m-2)^2;
solmeth = 2;
restart = 0;
prob = 1;
conv = 10^-5;
para = 4; %%%%% ARG %%%%%%

%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);hs =X(2)-X(1);
T = linspace(0,1,k);ht = T(2)-T(1);

vec = helpvector(m);

disk = ht^2/(hs^2);


%%%%%%%%%% TODO %%%%%%%%%%
%%% Endre alle funksjoner til å ta ut dimensjonsdata automatisk der
%%% mulig%%%
%%% Skrive bølge, maxwell og varme sammen så mye som mulig %%%
%%% Legge til beskrivelse til ferdige funksjoner %%%
%%% 



[ U0,V0,F1,F2,G1,G2,correctsolution] = getWaveTestFunctions( prob,X,T );

%%%% Matriser og vectorer %%%%
A = 1/hs^2*gallery('poisson', m-2);
%Atilde = [A,sparse((m-2)^2,(m-2)^2);sparse((m-2)^2,(m-2)^2),speye((m-2)^2)]; %% %%% IKKE BRUK!!!!%%


Jtilde = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)];
Atilde = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-A,sparse((m-2)^2,(m-2)^2)]; %Test

%Jtilde = [speye((m-2)^2),sparse((m-2)^2,(m-2)^2);sparse((m-2)^2,(m-2)^2),-speye((m-2)^2)];

U0tilde = [V0(vec);U0(vec)];
U0tilde1 = [U0(vec);V0(vec)];

F1tilde = [sparse((m-2)^2,1);F1(vec)];
F2tilde = [sparse((m-2)^2,1);F2(vec)];


%%%%TEST%%%%
%A = 1/hs^2*gallery('poisson', m-2);
%A = Atilde;
%J = Jtilde;
%v = U0tilde;

%[V,H,hn]=Arnoldi(A,v,4,conv);

%[Vtilde,Htilde] = SymplecticblockLanczos(A,J, V(:,1:4) , (m-2)^2 );




%%%%TEST%%%% Fungerer!
% Utemp = zeros(2*(m-2)^2,k);
% for i = 1:k
%     Utemp(:,i) = expm(Atilde*T(i))*U0tilde;
% end
%e1 = zeros(2*(m-2)^2,1); e1(1) = 1;
%eend = zeros(2*(m-2)^2,1); eend(end) = 1;
%allss = ones(2*(m-2)^2,1)/(2*(m-2)^2);
%q = [U0tilde,F1tilde,F2tilde,-Jtilde*U0tilde,-Jtilde*F1tilde,-Jtilde*F2tilde]; % HVa er q???
% Q = [U0tilde,F1tilde,F2tilde]; % HVa er q???
% Qtilde = mbsgs(Q,Jtilde);
% %q'*Jtilde*q == [0,1;-1,0]
% V = [Qtilde(:,1:1),-J*Qtilde(:,1:1)];
% [Vtilde,Htilde] = SymplecticblockLanczos(Atilde,Jtilde, V, n );
% %Vtilde = Vtilde/norm(U0tilde,2);
% Utemp = zeros(2*(m-2)^2,k);
% for i = 1:k
%     Utemp(:,i) = sum(Vtilde*expm(Htilde*T(i))*Vtilde'*Jtilde*V,2);
% end

%[Htilde,Vtilde,v] = SymplecticLanczosMethod(Atilde,Jtilde,U0tilde,n);
V = U0tilde1; [Vm,Hm,hm] = Arnoldi(Atilde,V,n,conv); Vm = Vm(:,1:n); V = Vm'*V;

%V = [U0tilde,-Jtilde*U0tilde]; [Vm,Hm] = SymplecticblockLanczos(Atilde,Jtilde, V , n/2 );V = Jtilde*V;
%V = U0tilde1; [Vm,Hm,~,~] = SymplecticLanczosMethod(Atilde,Jtilde,V,n/2);  V = Vm'*U0tilde1;
%V = [sparse(n/2,n/2),-speye(n/2);speye(n/2),sparse(n/2,n/2)]*Vm'*Jtilde*U0tilde;

%V = [sparse(n/2,n/2),-speye(n/2);speye(n/2),sparse(n/2,n/2)]*Vm'*Jtilde*U0tilde;%V = -[sparse((m-2)^2,(m-2)^2),-speye((m-2)^2);speye((m-2)^2),sparse((m-2)^2,(m-2)^2)]*Vm'*U0tilde;
%[Vm, Hm] = symplecticramSchmidt( Atilde,Jtilde, n/2 ); V = [sparse((m-2)^2,(m-2)^2),-speye((m-2)^2);speye((m-2)^2),sparse((m-2)^2,(m-2)^2)]*Vm'*Jtilde*U0tilde1;
%norm(U0tilde,2)
%%%%%%%%%%%%%%%%%%%%%%[Utemp,~] = KPM(Jtilde,Atilde,Jtilde*Atilde*U0tilde,1,k,n,2*(m-2)^2,ht,conv,restart);
Utemp = zeros(2*(m-2)^2,k);
for i = 1:k
    Utemp(:,i) = Vm*expm(Hm*T(i))*V;
end
%max(max(Utemp))


% Utemp = zeros(2*(m-2)^2,k);
% for i = 1:k
%     Utemp(:,i) = sum(Vm*expm(Hm*T(i))*Vm'*V,2);
% end
% max(max(Utemp))

% Utemp1 = zeros(2*(m-2)^2,k);
% for i = 1:k
%     Utemp1(:,i) = expm(Atilde*T(i))*U0tilde1;
% end
% max(max(Utemp1))
%max(Utemp1(:,1))/max(Utemp(:,1))
%Utemp = Utemp*max(max(Utemp1))/max(max(Utemp));


U = zeros(m^2,k);
U(vec,:) = Utemp(1:(m-2)^2,:);
% 
% U1 = zeros(m^2,k);
% U1(vec,:) = Utemp1(1:(m-2)^2,:);

%video(U,m,k,0.05)
%video(U1,m,k,0.05)
%video(U-U1,m,k,0.05)
%max(max(Utemp))
%max(max(correctsolution))
%Utemp = integrate(H,F,n,k,ht)



% if solmeth == 1
%     tic;
%     [Utemp1,iter1] = KPM(Jtilde,Atilde,Jtilde*Atilde*U0tilde,1,k,n,2*(m-2)^2,ht,conv,restart);
%     [Utemp2,iter2] = KPM(Jtilde,Atilde,F1tilde,G1,k,n,2*(m-2)^2,ht,conv,restart);
%     [Utemp3,iter3] = KPM(Jtilde,Atilde,F2tilde,G2,k,n,2*(m-2)^2,ht,conv,restart);
%     Utemp = Utemp1 + Utemp2 + Utemp3;
%     utdata(1) = iter1+iter2+iter3;
%     utdata(2) = toc;
% elseif solmeth == 2
%     utdata(1) = 0;
%     tic;
%     Utemp = integrate(Jtilde*Atilde,Jtilde*Atilde*U0tilde*ones(1,k)+F1tilde*G1+F2tilde*G2,2*(m-2)^2,k,ht);
%     utdata(2) = toc;
% end

%U = zeros(m^2,k);
%U(vec,:) = Utemp(1:(m-2)^2,:);
%U(vec,:) = U(vec,:) + U0(vec)*ones(1,k);

%V = zeros(m^2,k);
%V(vec,:) = Utemp((m-2)^2+1:end,:);
%V(vec,:) = V(vec,:) + U0(vec)*ones(1,k);

if 0
    video(U,m,k,0.05)
    %video(V,m,k,0.05)
    %video(correctsolution,m,k,0.05)
    %video(U-correctsolution,m,k,0.05)
    energy(Jtilde*Atilde,Utemp);
end

%energy(Jtilde*Atilde,Utemp-correctsolution(vec,(m)^2+vec,:)); % Fungerer
%ikke fordi den deriverte ikke blir tatt med i beregningene av energi


utdata(3) = max(max(abs(U-correctsolution)));
utdata
%utdata(3)
toc;