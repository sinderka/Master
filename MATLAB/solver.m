function [utdata] = solver(m,k,n,eqn,prob,solmeth,conv,para)

% m:     number of steps in each spacial direction (20)
% k:     number of steps in time (20)
% n:     restartvariable (20)
% eqn:   which equation to solve: heat, wave, maxwell (heat)
% prob:  test equation to solve (1)
% solmet:which solution method to use (restarted Krylov)
% conv:  Convergence criterion (10^-5)
% para:  number of processors used(4)

%%% default values
if nargin == 0
    m = 20;
    k = 20;
    n = 20;
    eqn = 'heat';
    prob = 3;
    solmeth = 3;
    conv = 10^-15;
    para = 4;
    
    
end

if strcmp(eqn,'heat')
    [utdata] = heatsolver(m,n,k,prob,solmeth,conv,para);
elseif strcmp(eqn,'wave')
    
elseif strcmp(eqn,'maxwell')
    
end

