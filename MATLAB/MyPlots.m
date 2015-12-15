function MyPlots
%%% All plots used in my master thesis
%%% Uncomment the rights sections to obtain pictures


%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI 
%plottool(m ,n ,k ,eqn         ,alg , int       ,restart, prob,conv,  para,data,type ,help,   name     ,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Første bildene på restart symp lanczos method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool(20,6,20,'semirandom',[-2,1,2],1,[-1,0,1]  ,1, 1e-14, 1    , 4 ,'table',0   ,'compareEnergy',1)
% plottool(20,6,20,'semirandom',[-2,1,2],1,[-1,0,1]  ,1, 1e-14, 1  ,   1 ,'table',    0   ,'compareIter',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Varying energy bildene på restart symp lanczos method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool(20,6,20,'semirandom',[-2,1,2],1,[-1,0,1]  ,2, 1e-14, 1    , 6 ,'table',0   ,'compareEnergy2',1)
% plottool(20,6,20,'semirandom',[-2,1,2],1,[-1,0,1]  ,2, 1e-14, 1  ,   1 ,'table',    0   ,'compareIter2',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Konstant energi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%energiTemp
% Andre bildene bildene på restart symp lanczos method
% energyTest(m,n,k,eqn,alg,integrator,restart,prob,conv,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver(20,6,20,'semirandom',2,1,0,1,1e-14,11) 
% figure(2) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
% saveit('energytestrestart0','t', 'energy2') 
% figure(5) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
% saveit('errortestrestart0','t', 'error2') 
% figure(7) 
% title(char({'semirandom', 'm=20','n=6', 'k=20', 'DI','trapezoidal rule','SLM'}))
% saveit('semirandomenergy1','t', 'energy') 
% 
% solver(20,6,20,'semirandom',2,1,1,1,1e-14,11) 
% figure(2) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','SLM'}))
% saveit('energytestrestart1','t', 'energy2') 
% figure(5) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','SLM'}))
% saveit('errortestrestart1','t', 'error2') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%energiTemp
% Tredje bildene bildene på restart symp lanczos method
% energyTest(m,n,k,eqn,~,integrator,restart,prob,conv,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver(20,6,20,'semirandom',1,1,0,1,1e-14,1) 
% figure(2) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('energyarnrestart0','t', 'energy2') 
% figure(5) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('errorarnrestart0','t', 'error2') 
% 
% solver(20,6,20,'semirandom',1,1,1,1,1e-14,1) 
% figure(2) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('energyarnrestart1','t', 'energy2') 
% figure(5) 
% title(char({'semirandom', 'm=20','n=4', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('errorarnrestart1','t', 'error2') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Non constant energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%energiTemp
% Andre bildene bildene på restart symp lanczos method
% energyTest(m,n,k,eqn,alg,integrator,restart,prob,conv,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver(20,6,20,'semirandom',2,1,0,2,1e-14,11) 
% figure(7) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
% saveit('energytestrestart02','t', 'energy2') 
% figure(5) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','SLM'}))
% saveit('errortestrestart02','t', 'error2') 
% figure(2) 
% title(char({'semirandom', 'm=20','n=6', 'k=20', 'DI','trapezoidal rule','SLM'}))
% saveit('semirandomenergy2','t', 'energy') 
% 
% solver(20,6,20,'semirandom',2,1,1,2,1e-14,11) 
% figure(7) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','SLM'}))
% saveit('energytestrestart12','t', 'energy2') 
% figure(5) 
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule'}))
% saveit('errortestrestart12','t', 'error2') 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%energiTemp
% Tredje bildene bildene på restart symp lanczos method
% energyTest(m,n,k,eqn,~,integrator,restart,prob,conv,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver(20,6,20,'semirandom',1,1,0,2,1e-14,1)
% figure(7)
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('energyarnrestart02','t', 'energy2')
% figure(5)
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('errorarnrestart02','t', 'error2')
% 
% solver(20,6,20,'semirandom',1,1,1,2,1e-14,1)
% figure(7)
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('energyarnrestart12','t', 'energy2')
% figure(5)
% title(char({'semirandom', 'm=20','n=6', 'k=20','restart=1', '\epsilon=1e-14','trapezoidal rule','KPM'}))
% saveit('errorarnrestart12','t', 'error2')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solver
% Figurer som er på andre rekke kappittel om integrasjon
% Viser hvordan energien ser ut til hver enkelt ser ut over tid. 
% Husk å kommentere inn plottet i funksjonen energy
%solver(m,n,k,eqn,alg,integrator,restart,prob,conv,para)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
solver(20,4,20,'wave',2,1,0,1,1e-14,1)
figure(2)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('energyovertimetrapezoidal','simulated time', 'energy1')
figure(11)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','trapezoidal rule'}))
saveit('errorovertimetrapezoidal','simulated time', 'error1')

solver(20,4,20^2,'wave',2,2,0,1,1e-14,1)
figure(2)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('energyovertimeeuler','t', 'energy1')
figure(11)
title(char({'wave', 'm=20^2','n=4', 'k=20','restart=0', '\epsilon=1e-14','forward Euler'}))
saveit('errorovertimeeuler','t', 'error1')

solver(20,4,20,'wave',2,3,0,1,1e-14,1)
figure(2)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('energyovertimemidpoint','t', 'energy1')
figure(11)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', '\epsilon=1e-14','midpoint rule'}))
saveit('errorovertimemidpoint','t', 'error1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solver
% Figurer som er på første rekke i kappittel om integrasjon
% Viser hvordan erroren konvergerer for de forkjellige integrasjonsmetodene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,[-1,10,20,40,80],'wave',2,[-2,1],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvtrap',1);
% saveit('intconvtrap', 'm=k', 'error1')
% plottool([-1,10,20,40,80],8,[-1,10^2,20^2,40^2,80^2],'wave',2,[-2,2],1,1, 1e-14, 1    , 3 ,'loglog',[1,-1]   ,'intconveul',1);
% saveit('intconveul', 'm=k', 'error1')
% plottool([-1,10,20,40,80],8,[-1,10,20,40,80],'wave',2,[-2,3],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvmid',1);
% saveit('intconvmid', 'm=k', 'error1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'resulttimem',1);
% plottool(20,8,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'resulttimek',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 1 ,'loglog',0   ,'resultiter',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,0,1, 1e-14, 1    , 1 ,'loglog',0   ,'resultiterr',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 3 ,'loglog',0   ,'resulterror',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,0,1, 1e-14, 1    , 3 ,'loglog',0   ,'resulterrorr',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,1, 1e-14, 1    , 4 ,'loglog',0   ,'resultenergy',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,0,1, 1e-14, 1    , 4 ,'loglog',0   ,'resultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'sresulttimem',1);
% plottool(20,8,[-1,10,20,40,80],'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 2 ,'loglog',0   ,'sresulttimek',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 1 ,'loglog',0   ,'sresultiter',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,0,1, 1e-14, 1    , 1 ,'loglog',0   ,'sresultiterr',1);
%%% plottool([-1,10,20,40,80],[-1,10,20,40,80],20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 5 ,'loglog',0   ,'sresulterror',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 5 ,'loglog',0   ,'sresulterror',1);
%%% plottool([-1,10,20,40,80],[-1,10,20,40,80],20,'semirandom',[-2,1,2,3],1,0,1, 1e-14, 1    , 5 ,'loglog',0   ,'sresulterrorr',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,0,1, 1e-14, 1    , 5 ,'loglog',0   ,'sresulterrorr',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,1, 1e-14, 1    , 4 ,'loglog',0   ,'sresultenergy',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,0,1, 1e-14, 1    , 4 ,'loglog',0   ,'sresultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vresulttimem',1);
% plottool(20,8,[-1,10,20,40,80],'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vresulttimek',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 1 ,'loglog',0   ,'vresultiter',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,0,2, 1e-14, 1    , 1 ,'loglog',0   ,'vresultiterr',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 3 ,'loglog',0   ,'vresulterror',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,0,2, 1e-14, 1    , 3 ,'loglog',0   ,'vresulterrorr',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,1,2, 1e-14, 1    , 4 ,'loglog',0   ,'vresultenergy',1);
% plottool([-1,10,20,40,80],8,20,'wave',[-2,1,2,3],1,0,2, 1e-14, 1    , 4 ,'loglog',0   ,'vresultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottools
% Figurer i "some interesting results"
% Arnoldi og slm vs tid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vsresulttimem',1);
% plottool(20,8,[-1,10,20,40,80],'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 2 ,'loglog',0   ,'vsresulttimek',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 1 ,'loglog',0   ,'vsresultiter',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,0,2, 1e-14, 1    , 1 ,'loglog',0   ,'vsresultiterr',1);
%%% plottool([-1,10,20,40,80],[-1,10,20,40,80],20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 5 ,'loglog',0   ,'vsresulterror',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 1    , 5 ,'loglog',0   ,'vsresulterror',1);
%%% plottool([-1,10,20,40,80],[-1,10,20,40,80],20,'semirandom',[-2,1,2,3],1,0,2, 1e-14, 1    , 5 ,'loglog',0   ,'vsresulterrorr',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,0,2, 1e-14, 1    , 5 ,'loglog',0   ,'vsresulterrorr',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,1,2, 1e-14, 2    , 4 ,'loglog',0   ,'vsresultenergy',1);
% plottool([-1,10,20,40,80],8,20,'semirandom',[-2,1,2,3],1,0,2, 1e-14, 2    , 4 ,'loglog',0   ,'vsresultenergyr',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
function saveit(name,xlab,ylab)
%legend(char(leg));
xlabel(xlab)
ylabel(ylab)
h = set(findall(gcf,'-property','FontSize'), 'Fontsize',18);
set(h,'Location','Best');
drawnow
    location = strcat('/home/shomeb/s/sindreka/Master/MATLAB/fig/',name);
    saveas(gcf,location,'fig');
    saveas(gcf,location,'jpeg');
end







