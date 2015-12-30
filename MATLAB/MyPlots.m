function MyPlots
%%% All plots used in my master thesis
% uncomment functions to make figures


% figures under chapter Energypreservation for SLM, constant energy
%SLMconstantenergy

% figures under chapter Energypreservation for SLM, variyng energy
%SLMvariyngenergy

% figures uunder chapter time integration methods
% convergence plots
%timeintegrationconvergence

% figures under chapter time integration methods
% names like energyovertimemidpoint or errorchangeretimetrapezoidal
%timeintegration

% figures under chapter K versus k
%Kversusk

% figures under chapter "the perfect restart variable"
%restartvariable

% figures under chapter "integrating over loong time"
integratinglongtime


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'resulttimem',1);
% plottool(20,8,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'resulttimek',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 1 ,'loglog',0   ,'resultiter',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,0,1, 1e-14, 1    , 1 ,'loglog',0   ,'resultiterr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 3 ,'loglog',0   ,'resulterror',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,0,1, 1e-14, 1    , 3 ,'loglog',0   ,'resulterrorr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 4 ,'loglog',0   ,'resultenergy',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,0,1, 1e-14, 1    , 4 ,'loglog',0   ,'resultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'sresulttimem',1);
% plottool(20,8,[-1,10,20,40,80],'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'sresulttimek',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 1 ,'loglog',0   ,'sresultiter',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,0,1, 1e-14, 1    , 1 ,'loglog',0   ,'sresultiterr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 5 ,'loglog',0   ,'sresulterror',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,0,1, 1e-14, 1    , 5 ,'loglog',0   ,'sresulterrorr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 4 ,'loglog',0   ,'sresultenergy',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,0,1, 1e-14, 1    , 4 ,'loglog',0   ,'sresultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vresulttimem',1);
% plottool(20,8,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vresulttimek',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 1 ,'loglog',0   ,'vresultiter',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,0,2, 1e-14, 1    , 1 ,'loglog',0   ,'vresultiterr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 3 ,'loglog',0   ,'vresulterror',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,0,2, 1e-14, 1    , 3 ,'loglog',0   ,'vresulterrorr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 4 ,'loglog',0   ,'vresultenergy',1);
% plottool([-1,10,20,40,80],8,1,1,20,'wave',[-2,1,2,3],1,0,2, 1e-14, 1    , 4 ,'loglog',0   ,'vresultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vsresulttimem',1);
% plottool(20,8,1,1,[-1,10,20,40,80],'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vsresulttimek',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 1 ,'loglog',0   ,'vsresultiter',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,0,2, 1e-14, 1    , 1 ,'loglog',0   ,'vsresultiterr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 5 ,'loglog',0   ,'vsresulterror',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,0,2, 1e-14, 1    , 5 ,'loglog',0   ,'vsresulterrorr',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 2    , 4 ,'loglog',0   ,'vsresultenergy',1);
% plottool([-1,10,20,40,80],8,1,1,20,'semirandom',[-2,1,2,3],1,0,2, 1e-14, 2    , 4 ,'loglog',0   ,'vsresultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Bilder i restartvariabel for tid og energi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool(40,[-1,4,6,8,10,20,40],1,1,40,'semirandom',[-2,1,2],1,1,1, 1e-14, 1    , 2 ,'semilogy',0   ,'restarttimeprob1',1);
% plottool(40,[-1,4,6,8,10,20,40],1,1,40,'semirandom',[-2,1,2],1,1,2, 1e-14, 1    , 2 ,'semilogy',0   ,'restarttimeprob2',1);
% plottool(40,[-1,4,6,8,10,20,40],1,1,40,'semirandom',[-2,1,2],1,0,1, 1e-14, 1    , 6 ,'semilogy',0   ,'restartenergyprob1',1);
% plottool(40,[-1,4,6,8,10,20,40],1,1,40,'semirandom',[-2,1,2],1,0,2, 1e-14, 1    , 6 ,'semilogy',0   ,'restartenergyprob2',1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

% Hvert kapitel/delkapitel som hører litt sammen har samme funksjon,
% funksjonsnavnet er et nøkkelord

%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para)

function SLMconstantenergy
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e5], 1    , 6 ,'loglog',0   ,'compareEnergy',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e5], 1    , 5 ,'loglog',0   ,'compareError',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e5], 1  ,   1 ,'semilogx',    0   ,'compareIter',1)

format shortEng; format compact; format shortEng; format compact; solver(20,6,1,1,20,'semirandom',2,1,0,1,1e-14,11)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
saveit('energytestrestart0','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
saveit('errortestrestart0','t', 'error2')
figure(2)
title(char({'semirandom', 'm=20','n=6', 'k=20', 'DI','trapezoidal rule','SLM'})) 
saveit('semirandomenergy1','t', 'energy')

format shortEng; format compact; format shortEng; format compact; solver(20,6,1,1,20,'semirandom',2,1,1,1,1e-6,11)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule','SLM'}))
saveit('energytestrestart2','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule','SLM'}))
saveit('errortestrestart2','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',2,1,1,1,1e-14,11)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','SLM'}))
saveit('energytestrestart1','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','SLM'}))
saveit('errortestrestart1','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',1,1,0,1,1e-14,1) 
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('energyarnrestart0','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('errorarnrestart0','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',1,1,1,1,1e-6,1)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule','KPM'}))
saveit('energyarnrestart2','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=4', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule','KPM'}))
saveit('errorarnrestart2','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',1,1,1,1,1e-14,1)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('energyarnrestart1','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=4', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('errorarnrestart1','t', 'error2')
end

function SLMvariyngenergy
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e5], 1    , 6 ,'loglog',0   ,'compareEnergy2',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e5], 1    , 5 ,'loglog',0   ,'compareError2',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e5], 1  ,   1 ,'semilogx',    0   ,'compareIter2',1)

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',2,1,0,2,1e-14,11)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
saveit('energytestrestart02','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
saveit('errortestrestart02','t', 'error2')
figure(2)
title(char({'semirandom', 'm=20','n=6', 'k=20', 'DI','trapezoidal rule','SLM'}))
saveit('semirandomenergy2','t', 'energy')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',2,1,1,2,1e-6,11)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule','SLM'}))
saveit('energytestrestart22','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule'}))
saveit('errortestrestart22','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',2,1,1,2,1e-14,11)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','SLM'}))
saveit('energytestrestart12','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule'}))
saveit('errortestrestart12','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',1,1,0,2,1e-14,1)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('energyarnrestart02','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('errorarnrestart02','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',1,1,1,2,1e-6,1)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule','KPM'}))
saveit('energyarnrestart22','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-6','trapezoidal rule','KPM'}))
saveit('errorarnrestart22','t', 'error2')

format shortEng; format compact; solver(20,6,1,1,20,'semirandom',1,1,1,2,1e-14,1)
format short
figure(7)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('energyarnrestart12','t', 'energy2')
figure(5)
title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
saveit('errorarnrestart12','t', 'error2')
end

function timeintegrationconvergence
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,1],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvtrap',1);
saveit('intconvtrap', 'm=k', 'error1')
plottool([-1,10,20,40,80],8,1,1,[-1,10^2,20^2,40^2,80^2],'wave',2,[-2,2],1,1, 1e-14, 1    , 3 ,'loglog',[1,-1]   ,'intconveul',1);
saveit('intconveul', 'm=k', 'error1')
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,3],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvmid',1);
saveit('intconvmid', 'm=k', 'error1')
end

function timeintegration


format shortEng; format compact; solver(20,4,1,1,20,'wave',2,1,0,1,1e-14,1)
format short
figure(2)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('energyovertimetrapezoidal','simulated time', 'energy1')
figure(11)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('errorovertimetrapezoidal','simulated time', 'error1')

format shortEng; format compact; solver(20,4,1,1,20^2,'wave',2,2,0,1,1e-14,1)
format short
figure(2)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('energyovertimeeuler','t', 'energy1')
figure(11)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('errorovertimeeuler','t', 'error1')

format shortEng; format compact; solver(20,4,1,1,20,'wave',2,3,0,1,1e-14,1)
format short
figure(2)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('energyovertimemidpoint','t', 'energy1')
figure(11)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('errorovertimemidpoint','t', 'error1')

format shortEng; format compact; solver(20,4,1,1,20,'wave',2,1,0,2,1e-14,1)
format short
figure(7)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('energychangtimetrapezoidal','simulated time', 'energy1')
figure(5)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('errorchangtimetrapezoidal','simulated time', 'error1')

format shortEng; format compact; solver(20,4,1,1,20^2,'wave',2,2,0,2,1e-14,1)
format short
figure(7)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('energychangtimeeuler','t', 'energy1')
figure(5)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('errorchangtimeeuler','t', 'error1')

format shortEng; format compact; solver(20,4,1,1,20,'wave',2,3,0,2,1e-14,1)
format short
figure(7)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('energychangtimemidpoint','t', 'energy1')
figure(5)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('errorchangtimemidpoint','t', 'error1')

% format shortEng; format compact; solver(20,4,1,1,20,'wave',3,1,0,2,1e-14,1)
% format short
% figure(7)
% title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
% saveit('energychangtimetrapezoidalt','simulated time', 'energy1')
% figure(5)
% title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
% saveit('errorchangtimetrapezoidalt','simulated time', 'error1')
% 
% format shortEng; format compact; solver(20,4,1,1,20^2,'wave',3,2,0,2,1e-14,1)
% format short
% figure(7)
% title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
% saveit('energychangtimeeulert','t', 'energy1')
% figure(5)
% title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
% saveit('errorchangtimeeulert','t', 'error1')
% 
% format shortEng; format compact; solver(20,4,1,1,20,'wave',3,3,0,2,1e-14,1)
% format short
% figure(7)
% title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
% saveit('energychangtimemidpointt','t', 'energy1')
% figure(5)
% title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
% saveit('errorchangtimemidpointt','t', 'error1')

end

function Kversusk
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,1,1e-14,1,2,'loglog',0,'Kversusktime',1)
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,1,1e-14,1,3,'loglog',0,'Kversuskerror',1)
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,1,1e-14,1,4,'loglog',0,'Kversuskenergy',1)
end

function restartvariable
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
plottool([-2,20,40,60],[-1,4,6,10,20,40],1,1,20,'wave',1,1,1,1,1e-14,4,2,'plot',0,'restarttime',1)
plottool([-2,20,40,60],[-1,4,6,10,20,40],1,1,20,'wave',1,1,1,1,1e-14,4,1,'plot',0,'restartiter',1)
plottool([-2,20,40,60],[-1,4,6,10,20,40],1,1,20,'wave',1,1,1,1,1e-14,4,3,'plot',0,'restarterror',1)
plottool([-2,20,40,60],[-1,4,6,10,20,40],1,1,20,'wave',1,1,1,1,1e-14,4,4,'plot',0,'restartenergy',1)
end

function integratinglongtime
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
plottool(20,6,[-1,1,5,10,20,40],1,1000,'wave',[-2,1,2,3],1,0,1,1e-14,4,3,'loglog',0,'longtime10',1)
plottool(20,6,[-1,1,5,10,20,40],1,1000,'wave',[-2,1,2,3],1,1,1,1e-14,4,3,'loglog',0,'longtime11',1)

plottool(20,6,[-1,1,5,10,20,40],1,1000,'wave',[-2,1,2,3],2,0,1,1e-14,4,3,'loglog',0,'longtime20',1)
plottool(20,6,[-1,1,5,10,20,40],1,1000,'wave',[-2,1,2,3],2,1,1,1e-14,4,3,'loglog',0,'longtime21',1)

plottool(20,6,[-1,1,5,10,20,40],1,1000,'wave',[-2,1,2,3],3,0,1,1e-14,4,3,'loglog',0,'longtime30',1)
plottool(20,6,[-1,1,5,10,20,40],1,1000,'wave',[-2,1,2,3],3,1,1,1e-14,4,3,'loglog',0,'longtime31',1)

end

function saveit(name,xlab,ylab)
%legend(char(leg));
xlabel(xlab)
ylabel(ylab)
h = set(findall(gcf,'-property','FontSize'), 'Fontsize',18);
set(h,'Location','Best');
pause(0.5)
drawnow
pause(0.5)
location = strcat('/home/shomeb/s/sindreka/Master/MATLAB/fig/',name);
saveas(gcf,location,'fig');
saveas(gcf,location,'jpeg');
end







