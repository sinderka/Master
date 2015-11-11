function [A] = getMatrix( m , hs, eqn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%% TODO Maxwell må fikses!

%%%%% Denne funksjonen skal fungere for varme, bølge og maxwell
if strcmp(eqn,'heat')
    A = -1/hs^2*gallery('poisson', m-2);
    
elseif strcmp(eqn, 'wave')
A = 1/hs^2*gallery('poisson', m-2);
A = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-A,sparse((m-2)^2,(m-2)^2)];
elseif strcmp(eqn,'maxwell') % DEnne MÅ testes!!!!!!!!!!
    %%%% ddx
T = 1/(2*hs)*gallery('tridiag',-ones(m-1,1),zeros(m,1),ones(m-1,1));
ddx = spalloc(m^3,m^3,4*m^2); % allokerer sikket for feil antall elementer
for i = 1:m^2
    ddx(1+(i-1)*m:i*m,1+(i-1)*m:i*m) = T;
end
%%%% ddy
ddy = spalloc(m^3,m^3,4*m^2);
for i = 1:m^2
    ddy(m*(i-1)+1+1,2*m+i) = 1;
    ddy(m*(i-1)+1+1,i) = -1;
    ddy(m*(i-1)+1,m+i) = 1;
    ddy(m*i,m+i) = -1;
end
%%%% ddz
ddz = spalloc(m^3,m^3,4*m^2);
for i = 1:m^2
    ddz(m*(i-1)+1+1,2*m^2+i) = 1;
    ddz(m*(i-1)+1+1,i) = -1;
    ddz(m*(i-1)+1,m^2+i) = 1;
    ddz(m*i,m^2+i) = -1;
end

%%%% Making curl!
my = 1; epsilon = 1; sigma = zeros(3*m^3,3*m^3);
T = [sparse(m^3,m^3),-ddz,ddy;
     ddz,sparse(m^3,m^3),-ddx;
     -ddy,ddx,sparse(m^3,m^3)];

A = [sparse(3*m^3,3*m^3), 1/my*T;1/epsilon*T,sigma];
    
end
end

