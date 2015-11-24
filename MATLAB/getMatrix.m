function [A] = getMatrix( m , hs, eqn )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%% TODO Maxwell må fikses!

%%%%% Denne funksjonen skal fungere for varme, bølge og maxwell
if strcmp(eqn,'heat')
    A = -1/hs^2*gallery('poisson', m-2);
    
elseif strcmp(eqn, 'wave')
    A = 1/hs^2*gallery('poisson', m-2);
    A = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-A,sparse((m-2)^2,(m-2)^2)];
elseif strcmp(eqn,'maxwell1D')
    G = 1/(2*hs)*gallery('tridiag',-ones(m-2,1),zeros(m-1,1),ones(m-2,1));
    
    
    %cols = 1/(2*hs)*[-ones(m-2,1),sparse(m-2,1),ones(m-2,1)];
    %G = sparse(1:m-2,1:m-2,cols,m-2,m-1);
    %spdiags([sparse(m-2,1),sparse(m-2,1),-ones(m-2,1),sparse(m-2,1),ones()])
    
    
    G(1,:) = [];G(1,1) = -1/hs; G(end,end+1) = 1/hs;A = -[sparse(m-2,m-2),G;-G',sparse(m,m)];
    
    A(1,m-1) = 0.5*A(1,m-1); A(m-2,end) = 0.5*A(m-2,end);
    %A = [sparse(m-2,m-2),G;G,sparse(m,m)];
    %G(1,2) = -1/hs; G(end,end-1) = -1/hs;
    %A = [sparse(m-2,m-2),G;G,sparse(m-2,m-2)];
    %isequal([sparse(m-2,m-2),speye(m-2);-speye(m-2),sparse(m-2,m-2)]*A,([sparse(m-2,m-2),speye(m-2);-speye(m-2),sparse(m-2,m-2)]*A)')
    a = 2;
    % elseif strcmp(eqn,'maxwell2D')
    %     T = 1/(2*hs)*gallery('tridiag',-ones(m-3,1),zeros(m-2,1),ones(m-3,1));
    %         %%%% ddx
    %     %T = 1/(2*hs)*gallery('tridiag',-ones(m-1,1),zeros(m,1),ones(m-1,1));
    %     ddx = spalloc(m,m^3,4*m^2); % allokerer sikket for feil antall elementer
    %     for i = 1:m
    %         ddx(1+(i-1)*m:i*m,1+(i-1)*m:i*m) = T;
    %     end
    %     %%%% ddy
    %     ddy = spalloc(m^3,m^3,4*m^2);
    %     for i = 1:m^2
    %         ddy(m*(i-1)+1+1,2*m+i) = 1;
    %         ddy(m*(i-1)+1+1,i) = -1;
    %         ddy(m*(i-1)+1,m+i) = 1;
    %         ddy(m*i,m+i) = -1;
    %     end
    %         A = [sparse(m^3,m^3),-ddy,ddx;
    %         -ddy,sparse(m^3,m^3),-ddx;
    %         ddx,ddx,sparse(m^3,m^3)];
    %
    % elseif strcmp(eqn,'maxwell3D') % DEnne MÅ testes!!!!!!!!!!
    %     %%%% ddx
    %     T = 1/(2*hs)*gallery('tridiag',-ones(m-1,1),zeros(m,1),ones(m-1,1));
    %     ddx = spalloc(m^3,m^3,4*m^2); % allokerer sikket for feil antall elementer
    %     for i = 1:m^2
    %         ddx(1+(i-1)*m:i*m,1+(i-1)*m:i*m) = T;
    %     end
    %     %%%% ddy
    %     ddy = spalloc(m^3,m^3,4*m^2);
    %     for i = 1:m^2
    %         ddy(m*(i-1)+1+1,2*m+i) = 1;
    %         ddy(m*(i-1)+1+1,i) = -1;
    %         ddy(m*(i-1)+1,m+i) = 1;
    %         ddy(m*i,m+i) = -1;
    %     end
    %     %%%% ddz
    %     ddz = spalloc(m^3,m^3,4*m^2);
    %     for i = 1:m^2
    %         ddz(m*(i-1)+1+1,2*m^2+i) = 1;
    %         ddz(m*(i-1)+1+1,i) = -1;
    %         ddz(m*(i-1)+1,m^2+i) = 1;
    %         ddz(m*i,m^2+i) = -1;
    %     end
    %
    %     %%%% Making curl!
    %     my = 1; epsilon = 1; sigma = zeros(3*m^3,3*m^3);
    %     T = [sparse(m^3,m^3),-ddz,ddy;
    %         ddz,sparse(m^3,m^3),-ddx;
    %         -ddy,ddx,sparse(m^3,m^3)];
    %
    %     A = [sparse(3*m^3,3*m^3), 1/my*T;1/epsilon*T,sigma];
    % elseif strcmp(eqn,'random')
    %     A = rand(2*(m-2)^2);
    %     A = 0.5*[sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)]*(A+A'+m^2*speye(2*(m-2)^2)); % eventuelt minus dette, men det har ikke noe å si!
    % elseif strcmp(eqn,'semirandom')
    %     try
    %         load('semirandomA.mat','A');
    %     catch
    %         A = -1;
    %     end
    %     if size(A,1) ~= 2*(m-2)^2
    %         A = rand(2*(m-2)^2);
    %         A = 0.5*[sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)]*(A+A'+m^2*speye(2*(m-2)^2)); % eventuelt minus dette, men det har ikke noe å si!
    %         save('semirandomA.mat','A');
    %     end
end
end

