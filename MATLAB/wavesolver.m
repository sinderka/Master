%function [ utdata ] = wavesolver( m,n,k,prob,solmeth,conv,para )

m = 50;
k = 400;
n = m^2;
solmeth = 3;
prob = 1;
%%% Initsiell data
utdata = zeros(1,3);
X = linspace(0,1,m);
hs =X(2)-X(1);
T = linspace(0,1,k);
ht = T(2)-T(1);
A = -1/hs^2*gallery('poisson', m-2);

disk = ht^2/(hs^2);

[U,V,F,correctsolution] = getWaveTestFunctions( prob );


if solmeth == 1
    
elseif solmeth == 2

elseif solmeth == 3 % it works!
    utdata(1) = 1;
    tic;
    U = doubleintegrate( U,V,F,A,ht,m,k );
    utdata(2) = toc;
    utdata(3) = max(max(max(abs(U-correctsolution))));
    
end

if 0
    video(U,m,k,0.05);
    video(correctsolution,m,k,0.05);
    video(U-correctsolution,m,k,0.05);
end
