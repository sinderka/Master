function [v,height,lastrelevant] = helpvector(m,eqn)
% Returns a list of indecis corresponding to non-edge points in a 1 and 2 D
% system.
%input:
% m: number of points in eqch spacial direction
% eqn: Decides if it is a 1 D og a 2 D system.
%Returns:
% v: a vector og sice m-2 (1D) or 2*(m-2)^2 corresponding to non edge
% points.
% height: m or m^2, dependant on the eqnation.


if strcmp(eqn,'maxwell1D')
    height = m;
    v = 2:m-1;
    lastrelevant = m-2; % Sikkert feil!
    return
end
    
height = m^2;
v = zeros((m-2)^2,1);
for qq = 0:m-3
    v(qq*(m-2)+1:qq*(m-2)+m-2) = (qq+1)*m+2:m-1 +(1+ qq)*m;
end
lastrelevant = (m-2)^2;
end