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
restartvariable

% figures under chapter "integrating over loong time"
%integratinglongtime

end

% Hvert kapitel/delkapitel som hører litt sammen har samme funksjon,
% funksjonsnavnet er et nøkkelord

%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para)

function SLMconstantenergy
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 6 ,'loglog',0   ,'compareEnergy',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 5 ,'loglog',0   ,'compareError',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e5], 1  ,   1 ,'semilogx',    0   ,'compareIter',1)

plottool(20,6,1,1,20,'wave',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 4 ,'loglog',0   ,'compareEnergyw',1)
plottool(20,6,1,1,20,'wave',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 3 ,'loglog',0   ,'compareErrorw',1)
plottool(20,6,1,1,20,'wave',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1  ,   1 ,'semilogx',    0   ,'compareIterw',1)


end

function SLMvariyngenergy
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 6 ,'loglog',0   ,'compareEnergy2',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 5 ,'loglog',0   ,'compareError2',1)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1  ,   1 ,'semilogx',    0   ,'compareIter2',1)

plottool(20,6,1,1,20,'wave',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 4 ,'loglog',0   ,'compareEnergy2w',1)
plottool(20,6,1,1,20,'wave',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , 3 ,'loglog',0   ,'compareError2w',1)
plottool(20,6,1,1,20,'wave',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1  ,   1 ,'semilogx',    0   ,'compareIter2w',1)

end

function timeintegrationconvergence
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,1],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvtrap',1);
saveit('intconvtrap', 'm=k', 'error1')
plottool([-1,10,20,40,80],8,1,1,[-1,10^2,20^2,40^2,80^2],'wave',2,[-2,2],1,1, 1e-14, 1    , 3 ,'loglog',[1,-1]   ,'intconveul',1);
saveit('intconveul', 'm^2=k', 'error1')
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,3],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvmid',1);
saveit('intconvmid', 'm=k', 'error1')

plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,1],1,2, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvtrap2',1);
saveit('intconvtrap', 'm=k', 'error1')
plottool([-1,10,20,40,80],8,1,1,[-1,10^2,20^2,40^2,80^2],'wave',2,[-2,2],1,2, 1e-14, 1    , 3 ,'loglog',[1,-1]   ,'intconveul2',1);
saveit('intconveul', 'm^2=k', 'error1')
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,3],1,2, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvmid2',1);
saveit('intconvmid', 'm=k', 'error1')
end

function timeintegration
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para)

format shortEng; format compact; solver(20,4,1,1,20,'wave',2,1,1,1,1e-14,1,1)
format short
figure(2)
title(char({'wave', 'm=20','n=4', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule'}))
saveit('energyovertimetrapezoidal','simulated time', 'energy1')
figure(11)
title(char({'wave', 'm=20','n=4', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule'}))
saveit('errorovertimetrapezoidal','simulated time', 'error1')

format shortEng; format compact; solver(20,4,1,1,20^2,'wave',2,2,1,1,1e-14,1,1)
format short
figure(2)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('energyovertimeeuler','t', 'energy1')
figure(11)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('errorovertimeeuler','t', 'error1')

format shortEng; format compact; solver(20,4,1,1,20,'wave',2,3,1,1,1e-14,1,1)
format short
figure(2)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('energyovertimemidpoint','t', 'energy1')
figure(11)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('errorovertimemidpoint','t', 'error1')

format shortEng; format compact; solver(20,4,1,1,20,'wave',2,1,1,2,1e-14,1,1)
format short
figure(7)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('energychangtimetrapezoidal','simulated time', 'energy1')
figure(5)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('errorchangtimetrapezoidal','simulated time', 'error1')

format shortEng; format compact; solver(20,4,1,1,20^2,'wave',2,2,1,2,1e-14,1,1)
format short
figure(7)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('energychangtimeeuler','t', 'energy1')
figure(5)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('errorchangtimeeuler','t', 'error1')

format shortEng; format compact; solver(20,4,1,1,20,'wave',2,3,1,2,1e-14,1,1)
format short
figure(7)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('energychangtimemidpoint','t', 'energy1')
figure(5)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('errorchangtimemidpoint','t', 'error1')
end

function Kversusk
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,1,1e-14,1,2,'loglog',0,'Kversusktime',1)
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,1,1e-14,1,3,'loglog',0,'Kversuskerror',1)
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,1,1e-14,1,4,'loglog',0,'Kversuskenergy',1)

plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,2,1e-14,1,2,'loglog',0,'Kversusktime2',1)
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,2,1e-14,1,3,'loglog',0,'Kversuskerror2',1)
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,2,1e-14,1,4,'loglog',0,'Kversuskenergy2',1)
end

function restartvariable
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,1,1e-14,4,2,'loglog',0,'restarttime',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,1,1e-14,4,1,'loglog',0,'restartiter',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,1,1e-14,4,3,'semilogy',0,'restarterror',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,1,1e-14,4,4,'semilogy',0,'restartenergy',1)

plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,1,1e-14,4,2,'loglog',0,'restarttimeSLM',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,1,1e-14,4,1,'loglog',0,'restartiterSLM',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,1,1e-14,4,3,'semilogy',0,'restarterrorSLM',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,1,1e-14,4,4,'semilogy',0,'restartenergySLM',1)


plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,2,1e-14,4,2,'loglog',0,'restarttime2',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,2,1e-14,4,1,'loglog',0,'restartiter2',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,2,1e-14,4,3,'semilogy',0,'restarterror2',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,2,1e-14,4,4,'semilogy',0,'restartenergy2',1)

plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,2,1e-14,4,2,'loglog',0,'restarttime2SLM',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,2,1e-14,4,1,'loglog',0,'restartiter2SLM',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,2,1e-14,4,3,'semilogy',0,'restarterror2SLM',1)
plottool([-2,20,40,60,80,100],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,2,1e-14,4,4,'semilogy',0,'restartenergy2SLM',1)
end

function integratinglongtime
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,0,1,1e-14,4,4,'loglog',[1,1],'longtime10',1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,1,1,1e-14,4,4,'loglog',0,'longtime11',1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,0,1,1e-14,4,4,'loglog',[1,2],'longtime20',1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,1,1,1e-14,4,4,'loglog',0,'longtime21',1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,0,1,1e-14,4,4,'loglog',[1,1],'longtime30',1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,1,1,1e-14,4,4,'loglog',0,'longtime31',1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,0,2,1e-14,4,4,'loglog',[1,1],'longtime102',1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,1,2,1e-14,4,4,'loglog',0,'longtime112',1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,0,2,1e-14,4,4,'loglog',[1,2],'longtime202',1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,1,2,1e-14,4,4,'loglog',0,'longtime212',1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,0,2,1e-14,4,4,'loglog',[1,1],'longtime302',1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,1,2,1e-14,4,4,'loglog',0,'longtime312',1)

end

function saveit(name,xlab,ylab)
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







