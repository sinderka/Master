%function [ utdata ] = wavesolver( m,n,k,prob,solmeth,conv,para )
% Få til en naiv implementasjon av 2D wave equation
clear
close all
m = 111;
k = 200;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k);
ht = T(2)-T(1);
A = -1/hs^2*gallery('poisson', m-2);

ant = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Feilen er at noe med kanten er galt!!!!! %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INIT U
U = zeros(m^2,k);
for i = 1:m
    for j = 1:m
        U((i)+(j-1)*m,1) = sin(pi*X(i))*sin(pi*X(j));
    end
end
% INIT F
%F = pi^2*sin(pi*X)'*ones(1,k);
%F = (2-X.*(1-X))'*sin(T);%
%F = -2*ones(m,k);
F = zeros(m^2,k);
%F = (1+pi^2)*sin(pi*X)'*exp(-T);
% INIT V
V = zeros(m^2,1);
%V = (2-X.*(1-X))';
%V = -sin(pi*X)';

helpvector = zeros((m-2)^2,1);
for i = 0:m-3
    helpvector(i*(m-2)+1:i*(m-2)+m-2) = (i+1)*m+2:m-1 +(1+ i)*m;
end

%U(:,1) = 1;
%for j = 0:m-1

%end
if 1
    U(helpvector,2) = U(helpvector,1) - ht * V(helpvector) + 0.5*ht^2*(A*U(helpvector,1) + F(helpvector,1));
    for i = 3:k
        %for j = 0:m-1
        U(helpvector,i) = 2*U(helpvector,i-1) - U(helpvector,i-2) + ht^2*(A*U(helpvector,i-1) + F(helpvector,i-1));
        %end
        mesh(reshape(U(:,i),m,m))
        axis([0,m,0,m,-1,1])
        %set ( handles.axes2, 'Color', [1 0 1] )
        drawnow
        pause(0.05)
    end
end
%mesh(U)

%correctsolution = sin(pi*X)'*cos(pi*T);
%correctsolution = sin(pi*X)'*cos(pi*T);
sol = @(t,x,y) sin(pi*x)*sin(pi*y)*cos(sqrt(2)*pi*t);
correctsolution = zeros(m^2,k);
for j = 1:k
    for i = 1:m
        for l = 1:m
            correctsolution(l+(i-1)*m,j) = sol(T(j),X(i),X(l));
        end
    end
end

%figure(2)
%mesh(correctsolution)
if 1
    for i = 3:k
        mesh(reshape(correctsolution(:,i),m,m))
        axis([0,m,0,m,-1,1])
        %set ( handles.axes2, )
        drawnow
        pause(0.05)
    end
end
%figure(3)
%mesh(U-correctsolution)
error = max(max(max(abs(U-correctsolution))))

if 1
    for i = 3:k
        mesh(reshape(U(:,i),m,m)-reshape(correctsolution(:,i),m,m))
        axis([0,m,0,m,min(min(U-correctsolution)),max(max(U-correctsolution))])
        %set ( handles.axes2, )
        drawnow
        pause(0.05)
    end
end
