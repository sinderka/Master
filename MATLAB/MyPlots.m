function MyPlots
%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI 
%plottool(m ,n ,k ,eqn         ,alg , int       ,restart, prob,conv,  para,data,type ,help,   name     ,save)
%start med -1: punkte på x aksen
%start med -2: forskjellige grapher

%plottool(20,6,20,'semirandom',[-2,1,2,3],1,[-1,0,1]  ,1, 1e-14, 1    , 4 ,'semilogy',0   ,'compareEnergy',1)
%plottool(20,6,20,'semirandom',[-2,1,2,3],1,[-1,0,1]  ,1, 1e-14, 1  ,   1 ,'plot',    0   ,'compareIter',1)
%plottool(20,6,20,'semirandom',[-2,1,2,3],1,[-1,0,1]  ,1, 1e-14, 1    , 6 ,'plot',0   ,'comparedifference',1)


%plottool(20,8,[-1,10,20,40,100],'wave',2,[-2,1,2,3],0,1, 1e-14, 1    , 6 ,'loglog',0   ,'semitime3integrateenergy',1)

%plottool(20,8,[-1,10,20,40,100],'wave',2,[-2,1,3],0,1, 1e-14, 1    , 6 ,'loglog',0   ,'semitime2integrateenergy',1)

%plottool(20,8,[-1,10,20,40,100],'semirandom',2,[-2,1,2,3],0,1, 1e-14, 1    , 5 ,'loglog',0   ,'semitime3integrateerror',1)

%plottool(20,8,[-1,10,20,40,100],'wave',2,[-2,1,3],1,1, 1e-14, 1    , 5 ,'loglog',0   ,'semitime2integrateerror',1)

%plottool(20,8,[-1,10,20,40,100],'wave',1,[-2,1,3],1,1, 1e-14, 1    , 5 ,'loglog',0   ,'semitime2integrateerror',1)

% Sammenligne energi for de forskjellige integratorene og de forskjellige
% metodene, for forskjellige PDE-er
%plottool(20,8,20,'wave',[-1,1,2],[-2,1,2,3],0,1, 1e-14, 1    , 3 ,'loglog',0   ,'waveerror1algint',1)
%plottool(20,8,20,'semirandom',[-1,1,2],[-2,1,2,3],0,1, 1e-14, 1    , 3 ,'loglog',0   ,'semiwaveerror1algint',1)
% Sammenligne error for de forskjellige integratorene og de forskjellige
% metodene, for forskjellige PDE-er
%plottool(20,8,20,'wave',[-1,1,2],[-2,1,3],0,1, 1e-14, 1    , 5 ,'loglog',0   ,'waveerror2algint',1)
%plottool(20,8,20,'semirandom',[-1,1,2],[-2,1,3],0,1, 1e-14, 1    , 5 ,'loglog',0   ,'semiwaveerror2algint',1)

% Viser hvordan energien ser ut til hver enkelt ser ut over tid. 
% Husk å kommentere inn plottet i funksjonen energy
%solver(m,n,k,eqn,alg,integrator,restart,prob,conv,para)
solver(20,4,20,'wave',2,1,0,1,1e-14,1)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', 'convergence criterion=1e-14','trapezoidal rule'}))
saveit('energyovertimetrapezoidal','time', 'Energy1')
solver(20,4,20,'wave',2,2,0,1,1e-14,1)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', 'convergence criterion=1e-14','forward Euler'}))
saveit('energyovertimeeuler','time', 'Energy1')
solver(20,4,20,'wave',2,3,0,1,1e-14,1)
title(char({'wave', 'm=20','n=4', 'k=20','restart=0', 'convergence criterion=1e-14','midpoint rule'}))
saveit('energyovertimemidpoint','time', 'Energy1')
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
