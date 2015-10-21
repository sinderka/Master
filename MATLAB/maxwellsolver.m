%function [ utdata ] = wavesolver( m,n,k,prob,solmeth,conv,para )
clear
close all
m = 3;
k = 20;
n = 1;
solmeth = 3;
prob = 1;
conv = 10^-15;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k);
ht = T(2)-T(1);
A = -1/hs^2*gallery('poisson', m-2);

%%%%%%%%%%%%%%%%%%%%%%%%% TODO %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% SJEKK OM DENNE FUNGERER FOR FORSKJELLIGE M %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
 
disk = ht^2/(hs^2);

[U,V,F,correctsolution] = getWaveTestFunctions( prob,m,k,X,T );