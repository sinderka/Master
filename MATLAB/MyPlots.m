%data: 1 == iter, 2 == time, 3 == error, 4 == energy, 5 == Difference in error between KPm and DI, 6 == differnce in energy between KPM and DI 
%plottool(m ,n ,k ,eqn         ,alg       ,restart, prob,conv,  para,data,type ,help,   name     ,save)
plottool(20,6,20,'semirandom',[-2,1,2,3],[-1,0,1]  ,1, 1e-14, 1    , 4 ,'semilogy',0   ,'compareEnergy',1)
plottool(20,6,20,'semirandom',[-2,1,2,3],[-1,0,1]  ,1, 1e-14, 1  ,   1 ,'plot',    0   ,'compareIter',1)
plottool(20,6,20,'semirandom',[-2,1,2,3],[-1,0,1]  ,1, 1e-14, 1    , 6 ,'plot',0   ,'comparedifference',1)