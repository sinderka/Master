function [v,height] = helpvector(m,eqn)

if strcmp(eqn,'maxwell1D')
    height = m;
    v = 2:m-1;
    return
end
    
height = m^2;
v = zeros((m-2)^2,1);
for qq = 0:m-3
    v(qq*(m-2)+1:qq*(m-2)+m-2) = (qq+1)*m+2:m-1 +(1+ qq)*m;
end
end