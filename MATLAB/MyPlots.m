function MyPlots
%%% All plots used in my master thesis
% uncomment functions to make figures


% figures under chapter Energypreservation for SLM, constant energy
SLMconstantenergy; display('DONE: SLMconstantenergy')

% figures under chapter Energypreservation for SLM, variyng energy
SLMvariyngenergy; display('DONE: SLMvariyngenergy')

% figures uunder chapter time integration methods
% convergence plots
timeintegrationconvergence; display('DONE: timeintegrationconvergence')

% figures under chapter time integration methods
% names like energyovertimemidpoint or errorchangeretimetrapezoidal
timeintegration; display('DONE: timeintegration')

% figures under chapter K versus k
Kversusk; display('DONE: Kversusk')

% figures under chapter "the perfect restart variable"
restartvariable; display('DONE: restartvariable')

% figures under chapter "integrating over loong time"
integratinglongtime; display('DONE: integratinglongtime')

% Figures under chapter "run time comparison"
runcomparison; display('DONE: runcomparison')

% Figures with SLM and long time 
SLMtime

end

% Hvert kapitel/delkapitel som hører litt sammen har samme funksjon,
% funksjonsnavnet er et nøkkelord

%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para)

function SLMconstantenergy
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'6','5','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergy','compareError','compareIter'},1)

plottool(20,6,1,1,20,'wave',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'4','3','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergyw','compareErrorw','compareIterw'},1)


end

function SLMvariyngenergy
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'6','5','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergy2','compareError2','compareIter2'},1)

plottool(20,6,1,1,20,'wave',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'4','3','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergy2w','compareError2w','compareIter2w'},1)

end

function timeintegrationconvergence
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,1],1,1,1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'intconvtrap'},1); 
saveit('intconvtrap', 'm=k', 'er_1')
pause(0.5)
plottool([-1,10,20,40,80],8,1,1,[-1,10^2,20^2,40^2,80^2],'wave',2,[-2,2],1,1, 1e-14, 1    , {'3'} ,{'loglog'},[1,-1]   ,{'intconveul'},1);
saveit('intconveul', 'm^2=k', 'er_1')
pause(0.5)
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,3],1,1, 1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'intconvmid'},1);
saveit('intconvmid', 'm=k', 'er_1')
pause(0.5)
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,1],1,2, 1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'intconvtrap2'},1);
saveit('intconvtrap', 'm=k', 'er_1')
pause(0.5)
plottool([-1,10,20,40,80],8,1,1,[-1,10^2,20^2,40^2,80^2],'wave',2,[-2,2],1,2, 1e-14, 1    , {'3'} ,{'loglog'},[1,-1]   ,{'intconveul2'},1);
saveit('intconveul', 'm^2=k', 'er_1')
pause(0.5)
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',2,[-2,3],1,2, 1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'intconvmid2'},1);
saveit('intconvmid', 'm=k', 'er_1')
end

function timeintegration
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,1)

savesolver(20,4,1,1,20,'wave',2,1,1,1,1e-14,1,{'energyovertimetrapezoidal','errorovertimetrapezoidal'},{'2','11'},'trapezoidal rule')

savesolver(20,4,1,1,20^2,'wave',2,2,1,1,1e-14,1,{'energyovertimeeuler','errorovertimeeuler'},{'2','11'},'forward Euler')

savesolver(20,4,1,1,20,'wave',2,3,1,1,1e-14,1,{'energyovertimemidpoint','errorovertimemidpoint'},{'2','11'},'midpoint rule')

savesolver(20,4,1,1,20,'wave',2,1,1,2,1e-14,1,{'energychangtimetrapezoidal','errorchangtimetrapezoidal'},{'7','5'},'trapezoidal rule')

savesolver(20,4,1,1,20^2,'wave',2,2,1,2,1e-14,1,{'energychangtimeeuler','errorchangtimeeuler'},{'7','5'},'forward Euler')

savesolver(20,4,1,1,20,'wave',2,3,1,2,1e-14,1,{'energychangtimemidpoint','errorchangtimemidpoint'},{'7','5'},'midpoint rule')

end

function Kversusk
% plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,0,1,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime0','Kversuskerror0','Kversuskenergy0'},1)

plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,0,2,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime20','Kversuskerror20','Kversuskenergy20'},1)

plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,1,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime','Kversuskerror','Kversuskenergy'},1)

plottool(20,6,9,[-1,40,20,10,4,2,1],[-1,10,20,40,100,200,400],'wave',[-2,1,2,3],1,1,2,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime2','Kversuskerror2','Kversuskenergy2'},1)
end

function restartvariable

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,1,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttime','restartiter','restarterror','restartenergy'},1)

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,1,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttimeSLM','restartiterSLM','restarterrorSLM','restartenergySLM'},1)

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,2,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttime2','restartiter2','restarterror2','restartenergy2'},1)

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,2,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttime2SLM','restartiter2SLM','restarterror2SLM','restartenergy2SLM'},1)
end

function integratinglongtime
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,0,1,1e-14,4,{'4'},{'loglog'},[1,1],{'longtime10'},1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime11'},1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,0,1,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime20'},1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime21'},1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,0,1,1e-14,4,{'4'},{'loglog'},[1,1],{'longtime30'},1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime31'},1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,0,2,1e-14,4,{'4'},{'loglog'},[1,1],{'longtime102'},1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],1,1,2,1e-14,4,{'4'},{'loglog'},0,{'longtime112'},1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,0,2,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime202'},1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],2,1,2,1e-14,4,{'4'},{'loglog'},0,{'longtime212'},1)

plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,0,2,1e-14,4,{'4'},{'loglog'},[1,1],{'longtime302'},1)
plottool(20,6,[-1,1,5,10,20,40,80],1,1000,'wave',[-2,1,2,3],3,1,2,1e-14,4,{'4'},{'loglog'},0,{'longtime312'},1)

end

function runcomparison
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,1  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimem'},1)

plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,0  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimem0'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimek'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimek0'},1)


plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,1  ,2,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimem'},1)

plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,0  ,2,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimem0'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1  ,2,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimek'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0  ,2,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimek0'},1)

end

function SLMtime

data = zeros(4,8);

data(1,:) = energyTest(20,4,20,1,'wave',0,1,1e-14);
data(2,:) = energyTest(20,4,20,10,'wave',0,1,1e-14);
data(3,:) = energyTest(20,4,20,100,'wave',0,1,1e-14);
data(4,:) = energyTest(20,4,20,1000,'wave',0,1,1e-14);

plot([1,10,100,1000],data(:,7));
%legend('')
saveit('SLMtime0','T_s','en_3')


data = zeros(4,8);

data(1,:) = energyTest(20,4,20,1,'wave',1,1,1e-14);
data(2,:) = energyTest(20,4,20,10,'wave',1,1,1e-14);
data(3,:) = energyTest(20,4,20,100,'wave',1,1,1e-14);
data(4,:) = energyTest(20,4,20,1000,'wave',1,1,1e-14);

plot([1,10,100,1000],data(:,7));
%legend('')
saveit('SLMtime1','T_s','en_3')

end

function savesolver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,name,figs,leg)

format shortEng; format compact; solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,1)
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
    [ylab,xlab,~,additionalInfo] = getLabels(1,m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,data);
    legend(leg)
    title(additionalInfo)
    saveit(name(i),xlab, ylab)
end
close all

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









