function MyPlots
%%% All plots used in my master thesis
% uncomment functions to make figures


% figures under chapter Energypreservation for SLM, constant energy
%SLMconstantenergy; display('DONE: SLMconstantenergy')

% figures under chapter Energypreservation for SLM, variyng energy
%SLMvariyngenergy; display('DONE: SLMvariyngenergy')

% figures uunder chapter time integration methods
% convergence plots
%timeintegrationconvergence; display('DONE: timeintegrationconvergence')

% figures under chapter time integration methods
% names like energyovertimemidpoint or errorchangeretimetrapezoidal
%timeintegration; display('DONE: timeintegration')

% figures under chapter K versus k
%Kversusk; display('DONE: Kversusk')

% figures under chapter "the perfect restart variable"
%restartvariable; display('DONE: restartvariable')

% figures under chapter "integrating over loong time"
%integratinglongtime; display('DONE: integratinglongtime')

% Figures under chapter "run time comparison"
%runcomparison; display('DONE: runcomparison')

% Figures with SLM and long time
SLMtime; display('DONE: SLMtime')

% Figures with SLM and epsilon
%SLMepsilon; display('DONE: SLMepsilon')

end

% Hvert kapitel/delkapitel som hører litt sammen har samme funksjon,
% funksjonsnavnet er et nøkkelord

%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para)

function SLMconstantenergy
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'6','5','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergy','compareError','compareIter'},1)

plottool(20,6,1,1,20,'wave',[-2,1,2],1,1,1,[-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'4','3','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergyw','compareErrorw','compareIterw'},1)


end

function SLMvariyngenergy
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
plottool(20,6,1,1,20,'semirandom',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'6','5','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergy2','compareError2','compareIter2'},1)

plottool(20,6,1,1,20,'wave',[-2,1,2],1,1  ,2, [-1,1e-14,1e-10,1e-6,1e-2,1e4], 1    , {'4','3','1'} ,{'loglog','loglog','semilogx'},0   ,...
    {'compareEnergy2w','compareError2w','compareIter2w'},1)

end

function timeintegrationconvergence
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)

plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,1,1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'intconv11'},0);
pause(1)
saveit('intconv11', 'm=k', 'er_1')
pause(1)
plottool([-1,10,20,40,80],8,1,1,[-1,10^2,20^2,40^2,80^2],'wave',[-2,1,2,3],2,0,1,1e-14, 1    , {'3'} ,{'loglog'},[1,-1]   ,{'intconv12'},0);
pause(1)
saveit('intconv12', 'm^2=k', 'er_1')
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,0,1,1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'intconv13'},0);
pause(1)
saveit('intconv13', 'm=k', 'er_1')


plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,3, 1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'invconv21'},0);
pause(1)
saveit('intconv21', 'm=k', 'er_1')
pause(1)
plottool([-1,10,20,40,80],8,1,1,[-1,10^2,20^2,40^2,80^2],'wave',[-2,1,2,3],2,1,3, 1e-14, 1    , {'3'} ,{'loglog'},[1,-1]   ,{'invconv22'},0);
pause(1)
saveit('intconv22', 'm^2=k', 'er_1')
pause(1)
plottool([-1,10,20,40,80],8,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],3,1,3, 1e-14, 1    , {'3'} ,{'loglog'},[1,-2]   ,{'invconv23'},0);
pause(1)
saveit('intconv23', 'm=k', 'er_1')

end

function timeintegration
% solver(m,n,simtime,K,k,eqn,alg,integrator,restart,prob,conv,para,1)

savesolver(20,4,1,1,20,'wave',2,1,1,1,1e-14,1,{'energyovertimetrapezoidal','errorovertimetrapezoidal'},{'2','11'},'trapezoidal rule')

savesolver(20,4,1,1,20^2,'wave',2,2,1,1,1e-14,1,{'energyovertimeeuler','errorovertimeeuler'},{'2','11'},'forward Euler')

savesolver(20,4,1,1,20,'wave',2,3,1,1,1e-14,1,{'energyovertimemidpoint','errorovertimemidpoint'},{'2','11'},'midpoint rule')

savesolver(20,4,1,1,20,'wave',2,1,1,3,1e-14,1,{'energychangtimetrapezoidal','errorchangtimetrapezoidal'},{'2','11'},'trapezoidal rule')

savesolver(20,4,1,1,20^2,'wave',2,2,1,3,1e-14,1,{'energychangtimeeuler','errorchangtimeeuler'},{'2','11'},'forward Euler')

savesolver(20,4,1,1,20,'wave',2,3,1,3,1e-14,1,{'energychangtimemidpoint','errorchangtimemidpoint'},{'2','11'},'midpoint rule')

end

function Kversusk
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)

plottool(20,6,5,[-1,100,50,20,10,5,2,1],[-1,1,2,5,10,20,50,100],'wave',[-2,1,2,3],1,0,1,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime0','Kversuskerror0','Kversuskenergy0'},1)

plottool(20,6,5,[-1,100,50,20,10,5,2,1],[-1,1,2,5,10,20,50,100],'wave',[-2,1,2,3],1,1,1,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime','Kversuskerror','Kversuskenergy'},1)

plottool(20,6,5,[-1,100,50,20,10,5,2,1],[-1,1,2,5,10,20,50,100],'wave',[-2,1,2,3],1,0,3,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime20','Kversuskerror20','Kversuskenergy20'},1)

plottool(20,6,5,[-1,100,50,20,10,5,2,1],[-1,1,2,5,10,20,50,100],'wave',[-2,1,2,3],1,1,3,1e-14,1,{'2','3','4'},{'loglog','loglog','loglog'},0,...
    {'Kversusktime2','Kversuskerror2','Kversuskenergy2'},1)
end

function restartvariable
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,1,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttime','restartiter','restarterror','restartenergy'},1)

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,1,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttimeSLM','restartiterSLM','restarterrorSLM','restartenergySLM'},1)

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',1,1,1,3,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttime2','restartiter2','restarterror2','restartenergy2'},1)

plottool([-2,20,40,60,80],[-1,4,6,10,20,40,80,100,120],1,1,20,'wave',2,1,1,3,1e-14,4,{'2','1','3','4'},{'loglog','loglog','semilogy','semilogy'},0,...
    {'restarttime2SLM','restartiter2SLM','restarterror2SLM','restartenergy2SLM'},1)
end

function integratinglongtime
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,0,1,1e-14,4,{'4'},{'loglog'},[1,1,1e-10],{'longtime10'},1)
plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime11'},1)

plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,0,1,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime20'},1)
plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime21'},1)

plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,0,1,1e-14,4,{'4'},{'loglog'},[1,1,1e-10],{'longtime30'},1)
plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,1,1,1e-14,4,{'4'},{'loglog'},0,{'longtime31'},1)

plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,0,3,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime102'},1)
plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],1,1,3,1e-14,4,{'4'},{'loglog'},0,{'longtime112'},1)

plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,0,3,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime202'},1)
plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],2,1,3,1e-14,4,{'4'},{'loglog'},0,{'longtime212'},1)

plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,0,3,1e-14,4,{'4'},{'loglog'},[1,2],{'longtime302'},1)
plottool(20,6,[-1,1,5,10,20,40,80,160],1,[-1,[1,5,10,20,40,80,160]*1000],'wave',[-2,1,2,3],3,1,3,1e-14,4,{'4'},{'loglog'},0,{'longtime312'},1)

end

function runcomparison
%plottool(m,n,simtime,K,k,'eqn',alg,int,restart,prob,conv,para,{data},'type',[help],{name},save,option)
plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,1  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimem'},1)

plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,0  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimem0'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimek'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0  ,1,1e-14, 1    , {'2'} ,{'loglog'},0,{'ccomparetimek0'},1)


plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,1  ,3,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimem'},1)

plottool([-1,10,20,40,80],6,1,1,20,'wave',[-2,1,2,3],1,0  ,3,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimem0'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1  ,3,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimek'},1)

plottool(20,6,1,1,[-1,10,20,40,80],'wave',[-2,1,2,3],1,0  ,3,1e-14, 1    , {'2'} ,{'loglog'},0,{'vcomparetimek0'},1)

end

function SLMtime
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save,option)
plottool(40,2,[-1,1,2,4,8,12,16,20,30,40],1,[-1,[1,2,4,8,12,16,20,30,40]*20],'wave',2,1,[-2,0,1],1,1e-6,4,{'1','2','3','4','7','8'},{'semilogx','semilogx','loglog','loglog','loglog','loglog'},[1,1,1e-10],{'SLMtimeiter','SLMtimetime','SLMtimeer1','SLMtimeen1','SLMtimeen3','SLMtimeen4'},1,2,0)
plottool(40,2,[-1,1,2^2,4^2,8^2,12^2,16^2,20^2,30^2,40^2],1,[-1,[1,2^2,4^2,8^2,12^2,16^2,20^2,30^2,40^2]*20],'wave',2,2,[-2,0,1],1,1e-6,4,{'1','2','3','4','7','8'},{'semilogx','semilogx','loglog','loglog','loglog','loglog'},[1,1,1e-10],{'SLMtimeiterf','SLMtimetimef','SLMtimeer1f','SLMtimeen1f','SLMtimeen3f','SLMtimeen4f'},1,2,0)


end

function SLMepsilon
%plottool(m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save,option)
plottool(20,4,1,1,20,'wave',2,1,1,[-2,1],[-1,1e-14,1e-10,1e-6,1e-2,1e4],4,{'1','7','8'},{'loglog','loglog','loglog'},0,{'SLMiter','SLMepsilonen3','SLMepsilonen4'},1,2,0)


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
    pause(1)
    saveit(name(i),xlab, ylab)
end
hold off

end

function saveit(name,xlab,ylab)
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









