%function [ utdata ] = wavesolver( m,n,k,prob,solmeth,conv,para )
clear
close all
m = 11;
k = 100;
n = 1;
solmeth = 2;
prob = 1;
conv = 10^-5;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k);
ht = T(2)-T(1);
A = -1/hs^2*gallery('poisson', m-2);

disk = ht^2/(hs^2);

[U,V,F,correctsolution] = getWaveTestFunctions( prob,m,k,X,T );


if solmeth == 1
    %U = zeros(m^2,k);
    tic;
    v = zeros((m-2)^2,1);
    %for i = 1:(m-2)^2
    v(1) = 1;
    v(end) = 1;
    [U,iter] = KPMwave( U,A,V,F,v,k,m,ht,(m-2)^2,conv );
    %U = U +Utemp;
    %v(i) = 0;
    %end
    utdata(2) = toc;
    utdata(1) = iter;
elseif solmeth == 2
        %U = zeros(m^2,k);
    tic;
    v = zeros((m-2)^2,1);
    %for i = 1:(m-2)^2
    v(1) = 1;
    v(end) = 1;
    [U,iter] = KPMwave( U,A,V,F,v,k,m,ht,n,conv );
    %U = U +Utemp;
    %v(i) = 0;
    %end
    utdata(2) = toc;
    utdata(1) = iter;
elseif solmeth == 3 % it works!
    utdata(1) = 1;
    tic;
    U = doubleintegrate( U,V,F,A,ht,m,k );
    utdata(2) = toc;
    
    
end
utdata(3) = max(max(max(abs(U-correctsolution))));
if 1
    video(U,m,k,0.05);
    video(correctsolution,m,k,0.05);
    video(U-correctsolution,m,k,0.05);
end
utdata