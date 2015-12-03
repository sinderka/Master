function MyPlots
%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI 
%plottool(m ,n ,k ,eqn         ,alg , int       ,restart, prob,conv,  para,data,type ,help,   name     ,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Første bildene på restart symp lanczos method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plottool(20,6,20,'semirandom',[-2,1,2,3],1,[-1,0,1]  ,1, 1e-14, 1    , 4 ,'semilogy',0   ,'compareEnergy',1)
% plottool(20,6,20,'semirandom',[-2,1,2,3],1,[-1,0,1]  ,1, 1e-14, 1  ,   1 ,'plot',    0   ,'compareIter',1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%energiTemp
% Andre bildene bildene på restart symp lanczos method
% energyTest(m,n,k,eqn,~,integrator,restart,prob,conv,~)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
energyTest(40,4,40,'wave',2,1,0,1,1e-14,11)
title(char({'wave', 'm=40','n=4', 'k=40','restart=0', 'convergence criterion=1e-14','trapezoidal rule'}))
saveit('energytestrestart0','time', 'Energy1')

energyTest(40,4,40,'wave',2,1,1,1,1e-14,11)
title(char({'wave', 'm=40','n=4', 'k=40','restart=1', 'convergence criterion=1e-14','trapezoidal rule'}))
saveit('energytestrestart1','time', 'Energy1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solver
% Figurer som er på andre rekke kappittel om integrasjon
% Viser hvordan energien ser ut til hver enkelt ser ut over tid. 
% Husk å kommentere inn plottet i funksjonen energy
%solver(m,n,k,eqn,alg,integrator,restart,prob,conv,para)
%%%%%%%%%%%%%%DISSE MANGLER LEGENDE OG SVART%%%%%%%%%%%%%%FARGE!!!!!!!!%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solver(20,4,20,'wave',2,1,0,1,1e-14,1)
% title(char({'wave', 'm=20','n=4', 'k=20','restart=0', 'convergence criterion=1e-14','trapezoidal rule'}))
% saveit('energyovertimetrapezoidal','time', 'Energy1')
% solver(20,4,20^2,'wave',2,2,0,1,1e-14,1)
% title(char({'wave', 'm=20','n=4', 'k=20','restart=0', 'convergence criterion=1e-14','forward Euler'}))
% saveit('energyovertimeeuler','time', 'Energy1')
% solver(20,4,20,'wave',2,3,0,1,1e-14,1)
% title(char({'wave', 'm=20','n=4', 'k=20','restart=0', 'convergence criterion=1e-14','midpoint rule'}))
% saveit('energyovertimemidpoint','time', 'Energy1')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solver
% Figurer som er på første rekke i kappittel om integrasjon
% Viser hvordan erroren konvergerer for de forkjellige integrasjonsmetodene
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plottool([-1,10,20,40,80],8,[-1,10,20,40,80],'wave',2,[-2,1],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvtrap',1);
%saveit('intconvtrap', 'm=k', 'error')
%plottool([-1,10,20,40,80],8,[-1,10^2,20^2,40^2,80^2],'wave',2,[-2,2],1,1, 1e-14, 1    , 3 ,'loglog',[1,-1]   ,'intconveul',1);
%saveit('intconveul', 'm=k', 'error')
%plottool([-1,10,20,40,80],8,[-1,10,20,40,80],'wave',2,[-2,3],1,1, 1e-14, 1    , 3 ,'loglog',[1,-2]   ,'intconvmid',1);
%saveit('intconvmid', 'm=k', 'error')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
function saveit(name,xlab,ylab)
%legend(char(leg));
xlabel(xlab)
ylabel(ylab)
h = set(findall(gcf,'-property','FontSize'), 'Fontsize',18);
set(h,'Location','Best');
    location = strcat('/home/shomeb/s/sindreka/Master/MATLAB/fig/',name);
    saveas(gcf,location,'fig');
    saveas(gcf,location,'jpeg');
end
