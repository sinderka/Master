function MyPlots
%%% All plots used in my master thesis
% uncomment functions to make figures

experiment; display('DONE: experiment');

% % 1
% % figures uunder chapter time integration methods
% % convergence plots
%timeintegrationconvergence; display('DONE: timeintegrationconvergence')

% % 2
% % Figures under chapter SLM energy
%longtime; display('DONE: longtime'); 

% 3
% Figures under chapter SLM perserved energy
%SLMperserveedenergy; display('DONE: SLMperservedenergy');

% 4
%stepbystep; display('DONE: stepbystep');
% % 5
% % figures under chapter K versus k
%Kversusk; display('DONE: Kversusk')
% 
% % 6
% % figures under chapter Energypreservation for SLM, constant energy
%changeeps; display('DONE: changeeps')

% % 7
% % figures under chapter "the perfect restart variable"
%restartvariable; display('DONE: restartvariable')

% 
% % 9
% % Figures in idea chapter
%ideaexpm; display('DONE: ideaexpm')
% 
% % 10
% % Figures in matalb expm
%matlabexpm; display('DONE: matlabexpm') 

% 
% % 8
% % Figures under chapter "run time comparison"
%runcomparison; display('DONE: runcomparison')

%pictures

%varyingenergy

end

% Hvert kapitel/delkapitel som hører litt sammen har samme funksjon,
% funksjonsnavnet er et nøkkelord

%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para)

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



% 1
function timeintegrationconvergence
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)

% with restart
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,1000;0,0,0]   ,{'intconv11','intener11'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10^2,20^2,40^2,80^2],'wave',[-2,1,2,3],2,1,1,1e-10, 1    ,[3,4] ,{'loglog','loglog'},[1,-1,1000;0,0,0]   ,{'intconv12','intener12'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,1,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,1000;0,0,0]   ,{'intconv13','intener13'},1,1,0,1);

% without 
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,1000;0,0,0]   ,{'intconv11r','intener11r'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10^2,20^2,40^2,80^2],'wave',[-2,1,2,3],2,0,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-1,1000;0,0,0]   ,{'intconv12r','intener12r'},1,1,0,1);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,0,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,1000;0,0,0]   ,{'intconv13r','intener13r'},1,1,0,1);

%exact
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,1,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,100;0,0,0]   ,{'exactconvtraper','exactconvtrapene'},1,1,0,3);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,1,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,100;0,0,0]   ,{'exactconvmider','exactconvmidene'},1,1,0,2);

plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,0,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,100;0,0,0]   ,{'exactconvtraperr','exactconvtrapener'},1,1,0,3);
plottool([-1,10,20,40,80],2,1,1,[-1,10,20,40,80],'wave',[-2,1,2],3,0,1,1e-10, 1    , [3,4] ,{'loglog','loglog'},[1,-2,100;0,0,0]   ,{'exactconvmiderr','exactconvmidener'},1,1,0,2);
end

% 2
function longtime
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
simtime = [1,2,4,8,12,16,20,40,80,100]; 
para = 1;

%without restart
plottool(20,200 ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],1,0,1,1e-6,para,[5,4],{'loglog','loglog'},[1,1,1e-15;1,1,1e-13],{'longtime2err','longtime2ene'},1,1,0,1)

%with restart
plottool(20,20  ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,1,1,1e-6,para,[5,4,1],{'loglog','loglog','loglog'},[1,1,2e-15;1,1,2e-13;0,1,1e-2],{'longtime2rerr','longtime2rene','longtime2rite'},1,1,0,1)
end

% 3
function SLMperserveedenergy
%very long time
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)
% simtime = [1,2,4,8,12,16,20,40,80,120,160,200,400,800,1200];
% plottool(20,200,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],1,0,1,1e-6,1,[5,4],{'loglog','loglog'},0,{'vlongerr','vlongene'},1,0,0,1)
% plottool(20,20,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,1,1,1e-6,1,[5,4],{'loglog','loglog'},0,{'vlongerrr','vlongener'},1,0,0,1)
% plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,[5,4],{'loglog','loglog'},[1,1,2e-15;1,1,2e-13],{'vlongerrrK','vlongenerK'},1,0,0,1)

% H3 or H4
simtime = [1,2,4,8,12,16,20,40,80,100];
m = 20; n = [-2,2,64]; K = 1; alg = [-2,2]; int = 3; restart = 1; prob = 1; para = 1; save = 1; PMint = 1; conv = 1e-6; data = [-3,7,8,9]; help = 0; type = {'loglog','loglog','loglog'};
plottool(m,m,[-1,simtime],K,[-1,20*simtime],'semirandom',alg,int,restart,prob,conv,para,data,type,help,{'SLMpes',},save,0,0,PMint)

end


% 4
function stepbystep

%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save,option,EToption,PMint)
simtime = [1,2,4,8,12,16,20,40,80,120,160,200,400];
m = 20; n = 200; K = 1; k = 200; alg = [-2,1,2]; int = 1; restart = 0; prob = 1; conv = 1e-6; para = 1; data =  [-3,10,11]; help = 0; save = 1; PMint = 1;
plottool(m,n,[-1,simtime],K,[-1,20*simtime],'semirandom',[-2,1],int,restart,prob,conv,para,data,{'loglog'},help,{'energswA'},save,0,0,PMint)
plottool(m,n,[-1,simtime],K,[-1,20*simtime],'semirandom',[-2,2],int,restart,prob,conv,para,data,{'loglog'},help,{'energswS'},save,0,0,PMint)

%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
simtime = [1,2,4,8,12,16,20,40,80,100]; 
%plottool(20,6,10,1,200,'semirandom',[-2,1,2],int,restart,prob,[-1,1:1:11],para,[5,4],{'semilogy','semilogy'},0,{'cierr2','ciene2'},save,0,0,PMint)
%plottool(20,20,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2],int,restart,prob,20,para,[5,4],{'loglog','loglog'},[1,1,2e-15;0,0,0],{'cierr1','ciene1'},save,0,0,PMint)
end

% 5
function Kversusk
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
simtime = [1,2,4,8,12,16,20,40,80,100]; 
plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,[2,5,4],{'loglog','loglog','loglog'},[0,0,0;1,1,2e-15;1,1,1e-13],{'Kversusktime0','Kversuskerror0','Kversuskenergy0'},1)
plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],3,1,1,1e-6,1,[2,5,4],{'loglog','loglog','loglog'},[0,0,0;1,1,2e-15;1,1,1e-13],{'Kversusktime','Kversuskerror','Kversuskenergy'},1)
end

% 6
function changeeps

plottool(20,20,100,1,2000,'wave',[-2,1,2],3,1,1,[-1,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [4,5,1] ,{'loglog','loglog','loglog'},0,{'lcompareEnergy','lcompareError','lcompareIter'},1,0,0,1)
plottool(20,20,100,1,2000,'semirandom',[-2,1,2],3,1,1,[-1,1e-12,1e-10,1e-8,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [4,5,1] ,{'loglog','loglog','loglog'},0,{'lcompareEnergyw','lcompareErrorw','lcompareIterw'},1,0,0,1)
end

% 7
function restartvariable

%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
% with restart
plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',1,3,1,1,1e-6,1,5,{'loglog'},0,{'reserrA'},1,0,0,1)
plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',2,3,1,1,1e-6,1,5,{'loglog'},0,{'reserrS'},1,0,0,1)

% without restart
plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',1,1,0,1,1e-6,1,5,{'loglog'},0,{'nerrorwA'},1,0,0,1)
plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',2,1,0,1,1e-6,1,5,{'loglog'},0,{'nerrorwS'},1,0,0,1)

% Time domain
plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',1,1,0,1,1e-6,1,5,{'loglog'},0,{'terrorwA'},1,0,0,1)
plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',2,1,0,1,1e-6,1,5,{'loglog'},0,{'terrorwS'},1,0,0,1)

plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',1,1,1,1,1e-6,1,5,{'loglog'},0,{'terrorwAr'},1,0,0,1)
plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',2,1,1,1,1e-6,1,5,{'loglog'},0,{'terrorwSr'},1,0,0,1)

end

% 8
function runcomparison 
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
vary = [4,8,12,16,20,40,80,100];
% without restart
plottool([-1,vary],200,100,1,20*100,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timem'},1,0,0,1)
plottool(20,200,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timek'},1,0,0,1)

% with restart
plottool([-1,vary],20,100,1,20*100,'semirandom',[-2,1,2,3],1,1,1,1e-6,1,2,{'loglog'},0,{'timemr'},1,0,0,1)
plottool(20,20,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],1,1,1,1e-6,1,2,{'loglog'},0,{'timekr'},1,0,0,1)
plottool(20,[-1,2,4,8,16,20,40,80,120,160,180,200,400,600],100,1,20*100,'semirandom',[-2,1,2,3],3,1,1,1e-6,1,2,{'loglog'},0,{'timekr1'},1,0,0,1)

% Windowing
plottool([-1,vary],20,100,100,20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timemt'},1,0,0,1)
plottool(20,20,[-1,vary],[-1,vary],20,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timekt'},1,0,0,1)

% Expm
plottool([-1,vary],20,100,1,20*100,'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timeme'},1,0,0,3)
plottool(20,20,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],1,0,1,1e-6,1,2,{'loglog'},0,{'timeke'},1,0,0,3)


end

% 9
function ideaexpm 
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)

plottool(20,20,[-1,1,2,4,8,12,16,20,30,40,80,100],1,[-1,[1,2,4,8,12,16,20,30,40,80,100]*20],'wave',[-2,1,2,3],1,0,1,1e-6,1,[3,4,2,5],{'loglog','loglog','loglog','loglog'},0,{'ideaerr20','ideaener20','ideatime20','ideacom20'},1,1,0,3)
%plottool(20,20,[-1,1,2,4,8,12,16,20,30,40,80,100],1,[-1,[1,2,4,8,12,16,20,30,40,80,100]*20],'wave',[-2,1,2,3],1,1,1,1e-6,1,[3,4,2],{'loglog','loglog','loglog'},0,{'ideaerr40','ideaener40','ideatime40'},1,1,0,3)
end

% 10
function matlabexpm
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{'data'},{'type'},[help],{'name'},save,option,EToption,PMint)

plottool(20,20,[-1,1,2,4,8,12,16,20,30,40,80,100],1,[-1,[1,2,4,8,12,16,20,30,40,80,100]*20],'wave',1,1,0,1,0,1,[3,4],{'loglog','loglog','loglog'},0,{'expmAerr','expmAener','expmAtime'},1,1,0,[-2,1,2,3])
plottool(20,20,[-1,1,2,4,8,12,16,20,30,40,80,100],1,[-1,[1,2,4,8,12,16,20,30,40,80,100]*20],'wave',2,1,0,1,0,1,[3,4],{'loglog','loglog','loglog'},0,{'expmSerr','expmSener','expmStime'},1,1,0,[-2,1,2,3])

end

function pictures

m = 20; n = 20;

plottool(20,20,[-1,10,20,30,40,50,60,70,80,90,100],1,[-1,20*[10,20,30,40,50,60,70,80,90,100]],'wave',[-2,1],1,0,1,1e-4,1,3,{'plot'},0,{'nice'},1,0,0,1)


energyTest(m,n,100,2000,'wave',1,0,1,1e-4,1,0,0,1,0)
[ylab,xlab,~,~] = getLabels(1,m,n,100,1,2000,'wave',1,1,0,1,1e-4,1,3,1);
figure(111);
legend('KPM(20)')
saveit('butterfly',xlab, ylab)
%
end


function varyingenergy
% % %plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,[data],{'type'},[help],{'name'},save,option,EToption,PMint)
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,3,1e-6, 1    , [3,4] ,{'loglog','loglog'},[1,-2,10;0,0,0]   ,{'varconv11','varener11'},1,1,0,1);
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,1,3,1e-6, 1    , [3,4] ,{'loglog','loglog'},[1,-2,10;0,0,0]   ,{'varconv13','varener13'},1,1,0,1);
% 
% % without 
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0,3,1e-6, 1    , [3,4] ,{'loglog','loglog'},[1,-2,10;0,0,0]   ,{'varconv11r','varener11r'},1,1,0,1);
% plottool([-1,10,20,40,80],20,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,0,3,1e-6, 1    , [3,4] ,{'loglog','loglog'},[1,-2,10;0,0,0]   ,{'varconv13r','varener13r'},1,1,0,1);
% 
% % Change eps
%plottool(20,20,100,1,2000,'wave',[-2,1,2],3,1,3,[-1,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [6,5,1] ,{'loglog','loglog','loglog'},0,{'varyEnergy','varyError','varyIter'},1,0,0,1)
%plottool(20,20,100,1,2000,'semirandom',[-2,1,2],3,1,2,[-1,1e-8,1e-6,1e-4,1e-2,1e-1,1e1,1e2,1e4], 1    , [6,5,1] ,{'loglog','loglog','loglog'},0,{'varyEnergyw','varyErrorw','varyIterw'},1,0,0,1)

% % 
% % Long time
% simtime = [1,2,4,8,12,16,20,40,80,100]; 
% para = 1;
% 
% % % without restart
% plottool(20,200 ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,0,2,1e-6,para,[5,6],{'loglog','loglog'},[1,1,1e-15;0,1,1e-13],{'varytime2err','varytime2ene'},1,1,0,1)
% % % 
% % % with restart
% plottool(20,20  ,[-1,simtime],1,[-1,20*simtime],'semirandom',[-2,1,2,3],3,1,2,1e-6,para,[5,6,1],{'loglog','loglog','loglog'},[1,1,2e-15;0,1,2e-13;0,1,1e-2],{'varytime2rerr','varytime2rene','varytime2rite'},1,1,0,1)
% simtime = [1,2,4,8,12,16,20,40,80,100]; 
% plottool(20,20,[-1,simtime],[-1,simtime],20,'semirandom',[-2,1,2,3],3,1,2,1e-6,1,[2,5,6],{'loglog','loglog','loglog'},[0,0,0;1,1,2e-15;1,1,1e-13],{'lversusktime0','lversuskerror0','lversuskenergy0'},1)
% % 
% 
% % restart variable
% plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',1,3,1,2,1e-6,1,3,{'loglog'},0,{'lreserrA'},1,0,0,1)
% plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',2,3,1,2,1e-6,1,3,{'loglog'},0,{'lreserrS'},1,0,0,1)
% 
% % without restart
% plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',1,3,0,2,1e-6,1,3,{'loglog'},0,{'lnerrorwA'},1,0,0,1)
% plottool([-2,10,20,40,60,80,100],[-1,2,4,8,16,20,40,80,120,160],10,1,20*10,'semirandom',2,3,0,2,1e-6,1,3,{'loglog'},0,{'lnerrorwS'},1,0,0,1)
% 
% % Time domain
% plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',1,3,0,2,1e-6,1,3,{'loglog'},0,{'lterrorwA'},1,0,0,1)
% plottool(20,[-1,8,16,20,40,80,120,160,200,300],[-2,[10,20,40,60,80,100]],1,[-2,20*[10,20,40,60,80,100]],'semirandom',2,3,0,2,1e-6,1,3,{'loglog'},0,{'lterrorwS'},1,0,0,1)
% 
% 
% % lag en runtime sak her!!!
vary = [4,8,12,16,20,40,80,100];
% % without restart
plottool([-1,vary],200,1,1,20,'semirandom',[-2,1,2,3],3,0,1,1e-6,1,2,{'loglog'},0,{'ltimem'},1,0,0,1)
plottool(20,200,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],3,0,1,1e-6,1,2,{'loglog'},0,{'ltimek'},1,0,0,1)

% with restart
plottool([-1,vary],20,1,1,20,'semirandom',[-2,1,2,3],3,1,2,1e-6,1,2,{'loglog'},0,{'ltimemr'},1,0,0,1)
plottool(20,20,[-1,vary],1,[-1,20*vary],'semirandom',[-2,1,2,3],3,1,2,1e-6,1,2,{'loglog'},0,{'ltimekr'},1,0,0,1)

end




function savesolver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,PMint,name,figs,leg) % blir brukt

format shortEng; format compact; solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,1,PMint)
format short

for i = 1:length(figs)
    text = strcat(strcat('figure(',figs(i)),')');
    eval(char(text))
    if str2num(char(figs(i))) == 2
        data = 4;
    elseif str2num(char(figs(i))) == 5
        data = 5;
    elseif str2num(char(figs(i))) == 7
        data = 6;
    elseif str2num(char(figs(i))) == 11
        data = 3;
    end
    [ylab,xlab,~,additionalInfo] = getLabels(1,m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,data,PMint);
    legend(leg)
    title(additionalInfo)
    pause(1)
    saveit(name(i),xlab, ylab)
    pause(1)
end
hold off

end 

function saveit(name,xlab,ylab) %blir brukt
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

function notUsed 
%function SLMtime
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save,option)
%plottool(40,2,[-1,1,2,4,8,12,16,20,30,40],1,[-1,[1,2,4,8,12,16,20,30,40]*20],'wave',2,1,[-2,0,1],1,1e-6,4,{'1','2','3','4','7','8'},{'semilogx','semilogx','loglog','loglog','loglog','loglog'},[1,1,1e-10],{'SLMtimeiter','SLMtimetime','SLMtimeer1','SLMtimeen1','SLMtimeen3','SLMtimeen4'},1,2,0)
%plottool(40,2,[-1,1,2^2,4^2,8^2,12^2,16^2,20^2,30^2,40^2],1,[-1,[1,2^2,4^2,8^2,12^2,16^2,20^2,30^2,40^2]*20],'wave',2,2,[-2,0,1],1,1e-6,4,{'1','2','3','4','7','8'},{'semilogx','semilogx','loglog','loglog','loglog','loglog'},[1,1,1e-10],{'SLMtimeiterf','SLMtimetimef','SLMtimeer1f','SLMtimeen1f','SLMtimeen3f','SLMtimeen4f'},1,2,0)
%end

%function SLMepsilon
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save,option)
%plottool(20,4,1,1,20,'wave',2,1,1,[-2,1],[-1,1e-14,1e-10,1e-6,1e-2,1e4],4,{'1','7','8'},{'loglog','loglog','loglog'},0,{'SLMiter','SLMepsilonen3','SLMepsilonen4'},1,2,0)
%end

% function integratinglongtime %Fjerne denne?
% %plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,0,1,1e-14,4,{'4'},{'loglog'},[1,1,1e-10],{'longtime10'},1)
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime11'},1)
% 
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,0,1,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime20'},1)
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime21'},1)
% 
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,0,1,1e-14,4,{'4'},{'loglog'},[1,1,1e-10],{'longtime30'},1)
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime31'},1)
% 
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,0,3,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime102'},1)
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,1,3,1e-14,4,{'4'},{'loglog'},0,{'longtime112'},1)
% 
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,0,3,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime202'},1)
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,1,3,1e-14,4,{'4'},{'loglog'},0,{'longtime212'},1)
% 
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,0,3,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime302'},1)
% plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,1,3,1e-14,4,{'4'},{'loglog'},0,{'longtime312'},1)
% 
% end

end





