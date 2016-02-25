function MyPlots
%%% All plots used in my master thesis
% uncomment functions to make figures
% -1
%experiment; display('DONE: experiment');

% pictures
% 0
%pictures; display('DONE: pictures');

% Constant energy
% 1
% convergence plots
%timeintegrationconvergence; display('DONE: timeintegrationconvergence')

% convergence with iota and i_r
% 2
%changeeps; display('DONE: changeeps'); 

% restartvariable
% 3
%restartvariable; display('DONE: restartvariable');

% longtime
% 4
%longtime; display('DONE: longtime')

% SLMperserveedenergy
% 5
%SLMperserveedenergy; display('DONE: SLMperserveedenergy')

% ideaexpm
% 6
%ideaexpm; display('DONE: ideaexpm')

% runcomparison
% 7
%runcomparison; display('DONE: runcomparison')

% varyingenergy
% 8
varyingenergy; display('DONE: varyingenergy')

end

% Hvert kapitel/delkapitel som hører litt sammen har samme funksjon,
% funksjonsnavnet er et nøkkelord

%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para)

% -1
function experiment
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)
%[-1,1,2,4,8,12,16,20,40,80,120,160,200,400,800,1200,1600,2000,4000,8000,12000,16000,20000]
%plottool(20,2,[-1,1,2,4,8,12,16,20,40,80,120],1,20,'wave',[-2,1,2,3],3,1,1,1e-14,1,[1,3,4],{'loglog','loglog','loglog'},0,{'exp1','exp2','exp3'},1,0,0,1)
%plottool(20,2,[-1,1,2,4,8,12,16,20,40,80,120],1,[-1,20*[1,2,4,8,12,16,20,40,80,120]],'wave',[-2,1,2,3],3,1,1,1e-14,1,[1,3,4],{'loglog','loglog','loglog'},0,{'exp4','exp5','exp6'},1,0,0,1)

%simtime = [1,2,4,8,12,16,20,40,80]; 
%plottool(20,[-2,12,20,30],10,1,20*10,'semirandom',1,3,1,1,[-1,1e-10,1e-8,1e-6,1e-5,1e-4,1e-3,1e-2],1,[5,4,2],{'loglog','loglog','loglog'},0,{'exp1','exp2','exp3'},1,0,0,1)

%plottool(40,[-2,12,20,30],10,1,20*10,'semirandom',1,3,1,1,[-1,1e-10,1e-8,1e-6,1e-5,1e-4,1e-3,1e-2],1,[5,4,2],{'loglog','loglog','loglog'},0,{'exp4','exp5','exp6'},1,0,0,1)

%plottool(20,[-2,4,12,20,50,100],10,1,20*10,'semirandom',1,3,1,1,[-1,1e-13,1e-12,1e-11,1e-10,1e-6,1e-2],1,[5,4,2],{'loglog','loglog','loglog'},0,{'exp4','exp5','exp6'},1,0,0,1)
simtime = [1,2,4,8,12,16,20,40,80,100]; 
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)
%plottool(20,20,[-1,simtime],1,[-1,20*simtime],'wave',[-2,1,2,3],1,0,1,1e-14,1,[3,4],{'loglog','loglog'},[1,1,5e-1;0,0,0],{'exp1','exp2'},1,0,0,1)

%plottool(20,20,[-1,simtime],1,[-1,20*simtime],'wave',[-2,1,2,3],1,0,3,1e-14,1,[3,4],{'loglog','loglog'},[0,1,5e-1;0,0,0],{'exp3','exp4'},1,0,0,1)
plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',1,1,1,1,1e-6,1,5,{'loglog'},0,{'terrorwAr'},1,0,0,1)
plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',2,1,1,1,1e-6,1,5,{'loglog'},0,{'terrorwSr'},1,0,0,1)
end

% 0
function pictures

m = 20; n = 20;

plottool(20,20,[-1,10,20,30,40,50,60,70,80,90,100],1,[-1,20*[10,20,30,40,50,60,70,80,90,100]],'wave',[-2,1],1,0,1,1e-4,1,3,{'plot'},0,{'nice'},1,0,0,1)


energyTest(m,n,100,2000,'wave',1,0,1,1e-4,1,0,0,1,0)
[ylab,xlab,~,~] = getLabels(1,m,n,100,1,2000,'wave',1,1,0,1,1e-4,1,3,1);
figure(111);
legend('KPM(20)')
saveit('butterfly',xlab, ylab)
close all
end

% 1
function timeintegrationconvergence
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)

% with restart
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,100;0,0,0]   ,{'intconv11','intener11'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10^2,20^2,40^2,80^2],'wave',[-2,1,2,3],2,1,1,1e-10, 1    ,[3,4] ,{'loglog','loglog'},[1,-1,100;0,0,0]   ,{'intconv12','intener12'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,1,1,1e-10, 1    , 3 ,{'loglog'},[1,-2,100]   ,{'intconv13'},1,1,0,1);

% without 
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0,1,1e-10, 1    , 3 ,{'loglog'},[1,-2,100]   ,{'intconv11r'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10^2,20^2,40^2,80^2],'wave',[-2,1,2,3],2,0,1,1e-10, 1    , 3 ,{'loglog'},[1,-1,100]   ,{'intconv12r'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,0,1,1e-10, 1    , 3 ,{'loglog'},[1,-2,100]   ,{'intconv13r'},1,1,0,1);

%exact
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,1,1,1e-10, 1    , 3 ,{'loglog'},[1,-2,10]   ,{'exactconvtraper'},1,1,0,3);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,1,1,1e-10, 1    , 3 ,{'loglog'},[1,-2,10]   ,{'exactconvmider'},1,1,0,2);

plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,0,1,1e-10, 1    , 3 ,{'loglog'},[1,-2,10]   ,{'exactconvtraperr'},1,1,0,3);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,0,1,1e-10, 1    , 3 ,{'loglog'},[1,-2,10]   ,{'exactconvmiderr'},1,1,0,2);
end

% 2
function changeeps

plottool(20,20,100,1,2000,'wave',[-2,1,2],3,1,1,[-1,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [4,5,1] ,{'loglog','loglog','loglog'},0,{'lcompareEnergy','lcompareError','lcompareIter'},1,0,0,1)
plottool(20,20,100,1,2000,'semirandom',[-2,1,2],3,1,1,[-1,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [4,5,1] ,{'loglog','loglog','loglog'},0,{'lcompareEnergyw','lcompareErrorw','lcompareIterw'},1,0,0,1)

%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
%%%%% NB commen inn line 62 and line 82 in PM.
simtime = [1,2,4,8,12,16,20,40,80,100]; 
para = 1; int = 3; prob = 1; restart = 1; save = 1; PMint = 1;
%plottool(20,[-2,6,6,12,12],10,1,200,'semirandom',[-2,1,2,1,2],int,1,prob,[-1,1:1:20],para,[5,4],{'semilogy','semilogy'},0,{'cierr2','ciene2'},save,0,0,PMint)
%plottool(20,20,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2],int,restart,prob,20,para,[5,4],{'loglog','loglog'},[1,1,2e-15;0,0,0],{'cierr1','ciene1'},save,0,0,PMint)
end

% 3
function restartvariable

%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
% with restart
plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',1,3,1,1,1e-6,1,5,{'loglog'},0,{'reserrA'},1,0,'NorthEast',1)
plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',2,3,1,1,1e-6,1,5,{'loglog'},0,{'reserrS'},1,0,'NorthEast',1)

% without restart
%plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',1,1,0,1,1e-6,1,5,{'loglog'},0,{'nerrorwA'},1,0,'NorthEast',1)
%plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',2,1,0,1,1e-6,1,5,{'loglog'},0,{'nerrorwS'},1,0,'NorthEast',1)

% Time domain without restart
plottool(20,[-1,8,16,20,40,80,120,160,200,300,400,600],[-2,[10,20,60,100]],1,[-2,20*[10,20,60,100]],'semirandom',1,1,0,1,1e-6,1,5,{'loglog'},0,{'terrorwA'},1,0,'SouthWest',1)
plottool(20,[-1,8,16,20,40,80,120,160,200,300,400,600],[-2,[10,20,60,100]],1,[-2,20*[10,20,60,100]],'semirandom',2,1,0,1,1e-6,1,5,{'loglog'},0,{'terrorwS'},1,0,'SouthWest',1)
% Time domain with restart
plottool(20,[-1,8,16,20,40,80,120,160,200,300,400,600],[-2,[10,20,60,100]],1,[-2,20*[10,20,60,100]],'semirandom',1,1,1,1,1e-6,1,5,{'loglog'},0,{'terrorwAr'},1,0,'NorthEast',1)
plottool(20,[-1,8,16,20,40,80,120,160,200,300,400,600],[-2,[10,20,60,100]],1,[-2,20*[10,20,60,100]],'semirandom',2,1,1,1,1e-6,1,5,{'loglog'},0,{'terrorwSr'},1,0,'NorthEast',1)
end

% 4
function longtime
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
simtime = [1,2,4,8,12,16,20,40,80,100]; 
para = 1;

%without restart
plottool(20,200 ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],1,0,1,1e-6,para,[5,4],{'loglog','loglog'},[1,1,1e-15;0,1,1e-13],{'longtime2err','longtime2ene'},1,1,0,1)

%with restart
plottool(20,20  ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,1,1,1e-6,para,[5,4,1],{'loglog','loglog','loglog'},[1,1,2e-15;0,1,2e-13;0,1,1e-2],{'longtime2rerr','longtime2rene','longtime2rite'},1,1,0,1)

% windowing
simtime = [1,2,4,8,12,16,20,40,80,100]; 
plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,[2,5,4],{'loglog','loglog','loglog'},[0,0,0;1,1,2e-10;0,1,1e-13],{'Kversusktime0','Kversuskerror0','Kversuskenergy0'},1)
plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],3,1,1,1e-6,1,[2,5,4],{'loglog','loglog','loglog'},[0,0,0;1,1,2e-15;0,1,1e-13],{'Kversusktime','Kversuskerror','Kversuskenergy'},1)
end


% 5
function SLMperserveedenergy
%very long time
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)
simtime = [1,2,4,8,12,16,20,40,80,120,160,200,400,800,1200];
plottool(20,200,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],1,0,1,1e-6,1,[5,4],{'loglog','loglog'},0,{'vlongerr','vlongene'},1,0,0,1)
plottool(20,20,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,1,1,1e-6,1,[5,4],{'loglog','loglog'},0,{'vlongerrr','vlongener'},1,0,0,1)
plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,[5,4],{'loglog','loglog'},[1,1,2e-10;1,1,2e-13],{'vlongerrrK','vlongenerK'},1,0,0,1)

% Energy in transformations
simtime = [1,2,4,8,12,16,20,40,80,120,160,200,400];
m = 20; n = 200; K = 1; k = 200; alg = [-2,1,2]; int = 1; restart = 0; prob = 1; conv = 1e-6; para = 1; data =  [-3,10,11]; help = 0; save = 1; PMint = 1;
plottool(m,n,[-1,simtime],K,[-1,20*simtime],'semirandom',[-2,1],int,restart,prob,conv,para,data,{'loglog'},help,{'energswA'},save,0,0,PMint)
plottool(m,n,[-1,simtime],K,[-1,20*simtime],'semirandom',[-2,2],int,restart,prob,conv,para,data,{'loglog'},help,{'energswS'},save,0,0,PMint)

% Residual energy
% simtime = [1,2,4,8,12,16,20,40,80,100];
% m = 20; n = [-2,2,64]; K = 1; alg = [-2,2]; int = 3; restart = 1; prob = 1; para = 1; save = 1; PMint = 1; conv = 1e-6; data = [-3,7,8,9]; help = 0; type = {'loglog','loglog','loglog'};
% plottool(m,m,[-1,simtime],K,[-1,20*simtime],'semirandom',alg,int,restart,prob,conv,para,data,type,help,{'SLMpes',},save,0,0,PMint)

end

% 6
function ideaexpm 
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)

% idea
plottool(20,20,[-1,1,2,4,8,12,16,20,30,40,80,100],1,[-1,[1,2,4,8,12,16,20,30,40,80,100]*20],'wave',[-2,1,2,3,1,2],1,0,1,1e-6,1,[3,4],{'loglog','loglog'},[1,1,4e-1;0,0,0],{'ideaerr20','ideaener20'},1,1,0,[-2,3,3,1,1,1])

% Matlabs expm
simtime = [1,2,4,8,12,16,20,40,80,120,160,200,400,800,1200];
plottool(20,20,[-1,simtime],1,[-1,simtime*20],'wave',[-2,1,2,1,2  ],1,0,1,1e-6,1,[3,4],{'loglog','loglog'},0                ,{'expmAerr','expmAener' },1,1,0,[-2,3,3,2,2])
end

% 7
function runcomparison 
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
vary = [4,8,12,16,20,40,80,100];
% without restart
plottool([-1,vary],200,100,1,20*100,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timem'},1,0,0,1)
plottool(20,200,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timek'},1,0,0,1)
%plottool(20,200,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timek1'},1,0,0,1)
plottool(20,[-1,2,4,8,16,20,40,80,120,200,400,600],100,1,20*100,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timek1'},1,0,0,1)

% with restart
plottool([-1,vary],20,100,1,20*100,'semirandom',[-2,1,2,3],1,1,1,1e-6,1,2,{'loglog'},0,{'timemr'},1,0,0,1)
plottool(20,20,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],1,1,1,1e-6,1,2,{'loglog'},0,{'timekr'},1,0,0,1)
plottool(20,[-1,2,4,8,16,20,40,80,120,200,400,600],100,1,20*100,'semirandom',[-2,1,2,3],3,1,1,1e-6,1,2,{'loglog'},0,{'timekr1'},1,0,0,1)

% Windowing
plottool([-1,vary],20,100,100,20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timemt'},1,0,0,1)
plottool(20,20,[-1,vary],[-1,vary],20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timekt'},1,0,0,1)

plottool([-1,vary],20,100,100,20,'semirandom',[-2,1,2,3],1,1,1,1e-6,1,2,{'loglog'},0,{'timemtr'},1,0,0,1)
plottool(20,20,[-1,vary],[-1,vary],20,'semirandom',[-2,1,2,3],1,1,1,1e-6,1,2,{'loglog'},0,{'timektr'},1,0,0,1)

% Expm
plottool([-1,vary],20,100,1,20*100,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timeme'},1,0,0,3)
plottool(20,20,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timeke'},1,0,0,3)


end

% 8
function varyingenergy
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
%convergence with restart
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,3,1e-6, 1    , 3 ,{'loglog'},[1,-2,1;0,0,0]   ,{'varconv11'},1,1,0,1);
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,1,3,1e-6, 1    , 3 ,{'loglog'},[1,-2,1;0,0,0]   ,{'varconv13'},1,1,0,1);
% 
% % convergence without restart
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0,3,1e-6, 1    , 3 ,{'loglog'},[1,-2,1;0,0,0]   ,{'varconv11r'},1,1,0,1);
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,0,3,1e-6, 1    , 3 ,{'loglog'},[1,-2,1;0,0,0]   ,{'varconv13r'},1,1,0,1);
%  
% % % % Change eps
% plottool(20,20,100,1,2000,'wave',[-2,1,2],3,1,3,[-1,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [6,5,1] ,{'loglog','loglog','loglog'},0,{'varyEnergy','varyError','varyIter'},1,0,0,1)
% plottool(20,20,100,1,2000,'semirandom',[-2,1,2],3,1,2,[-1,1e-8,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [6,5,1] ,{'loglog','loglog','loglog'},0,{'varyEnergyw','varyErrorw','varyIterw'},1,0,0,1)

% restartvariable
% with restart
plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'wave',1,3,1,3,1e-6,1,5,{'loglog'},0,{'lreserrA'},1,0,'NorthEast',1)
plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'wave',2,3,1,3,1e-6,1,5,{'loglog'},0,{'lreserrS'},1,0,'NorthEast',1)

% % without restart
plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'wave',1,3,0,3,1e-6,1,5,{'loglog'},0,{'lnerrorwA'},1,0,'SouthWest',1)
plottool([-2,10,20,60,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'wave',2,3,0,3,1e-6,1,5,{'loglog'},0,{'lnerrorwS'},1,0,'SouthWest',1)

% Time domain
plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,60,100]],1,[-2,20*[10,20,60,100]],'semirandom',1,1,0,2,1e-6,1,5,{'loglog'},0,{'lterrorwA'},1,0,'SouthWest',1)
plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,60,100]],1,[-2,20*[10,20,60,100]],'semirandom',2,1,0,2,1e-6,1,5,{'loglog'},0,{'lterrorwS'},1,0,'SouthWest',1)




% % % Energy and error
% simtime = [1,2,4,8,12,16,20,40,80,100]; 
% para = 1;
% %without restart
% plottool(20,200 ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,0,2,1e-6,para,[5,6],{'loglog','loglog'},[1,1,1e-15;0,1,1e-13],{'vlongtime2err','vlongtime2ene'},1,1,0,1)
% %with restart
% plottool(20,20  ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,1,2,1e-6,para,[5,6,1],{'loglog','loglog','loglog'},[1,1,2e-15;0,1,2e-13;0,1,1e-2],{'vlongtime2rerr','vlongtime2rene','vlongtime2rite'},1,1,0,1)
% 
% % % Windowing
% simtime = [1,2,4,8,12,16,20,40,80,100]; 
% plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],3,0,2,1e-6,1,[5,6],{'loglog','loglog'},[0,0,0;0,1,2e-15;0,1,1e-13],{'lversuskerror0','lversuskenergy0'},1)
% plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],3,1,2,1e-6,1,[5,6],{'loglog','loglog'},[0,0,0;0,1,2e-15;0,1,1e-13],{'lversuskerror0r','lversuskenergy0r'},1)
% 
% % Long time
% simtime = [1,2,4,8,12,16,20,40,80,120,160,200,400,800,1200];
% plottool(20,200 ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,0,2,1e-6,para,[5,6],{'loglog','loglog'},[1,1,1e-15;0,1,1e-13],{'longererr','longerene'},1,1,0,1)
% plottool(20,20  ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,1,2,1e-6,para,[5,6,1],{'loglog','loglog','loglog'},[1,1,2e-15;0,1,2e-13;0,1,1e-2],{'longererrr','longerener','longeriter'},1,1,0,1)
% 
% % For wave
% simtime = [1,2,4,8,12,16,20,40,80,100]; 
% plottool(20,200,[-1,simtime],1,[-1,20*simtime],'wave',[-2,1,2,3],3,0,3,1e-6,1,[3,4],{'loglog','loglog'},0,{'vwaveerr','vwaveene'},1,0,0,1)
% plottool(20,20,[-1,simtime],1,[-1,20*simtime],'wave',[-2,1,2,3],3,1,3,1e-6,1,[3,4,1],{'loglog','loglog','loglog'},0,{'vwavererr','vwaverene','vwaveiter'},1,0,0,1)
% 
% % % runtime
% vary = [4,8,12,16,20,40,80,100];
% % without restart
% plottool([-1,vary],200,1,1,20,'semirandom',[-2,1,2,3],3,0,1,1e-6,1,2,{'loglog'},0,{'ltimem'},1,0,0,1)
% plottool(20,200,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],3,0,1,1e-6,1,2,{'loglog'},0,{'ltimek'},1,0,0,1)
% plottool(20,[-1,2,4,8,16,20,40,80,120,200,400,600],100,1,20*100,'semirandom',[-2,1,2,3],3,0,2,1e-6,1,2,{'loglog'},0,{'ltimek1'},1,0,0,1)
% 
% % with restart
% plottool([-1,vary],20,1,1,20,'semirandom',[-2,1,2,3],3,1,2,1e-6,1,2,{'loglog'},0,{'ltimemr'},1,0,0,1)
% plottool(20,20,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],3,1,2,1e-6,1,2,{'loglog'},0,{'ltimekr'},1,0,0,1)
% plottool(20,[-1,2,4,8,16,20,40,80,120,200,400,600],100,1,20*100,'semirandom',[-2,1,2,3],3,1,2,1e-6,1,2,{'loglog'},0,{'ltimekr1'},1,0,0,1)
end



function saveit(name,xlab,ylab)
xlabel(xlab)
ylabel(ylab)
h = set(findall(gcf,'-property','FontSize'), 'Fontsize',18);
set(h,'Location','Best');
pause(0.5)
drawnow
pause(0.5)
location = strcat('/home/shomeb/s/sindreka/Master/MATLAB/fig/',char(name));
saveas(gcf,location,'fig');
saveas(gcf,location,'jpeg');

end



