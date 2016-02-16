function [A] = getMatrix( m , hs, eqn )
% returns a matrix dependantn on what problems is to be sovled
%input:
% m: number of points in space
% hs: stempelngth in space
% eqn: a string coresponding to an equation
%returns:
% a matrix

%%%%% Denne funksjonen skal fungere for varme, bølge og maxwell
if strcmp(eqn,'heat')
    A = -1/hs^2*gallery('poisson', m-2);
    
elseif strcmp(eqn, 'wave')
    A = 1/hs^2*gallery('poisson', m-2);
    A = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-A,sparse((m-2)^2,(m-2)^2)];
elseif strcmp(eqn,'maxwell1D')
    G = 1/(2*hs)*gallery('tridiag',-ones(m-2,1),zeros(m-1,1),ones(m-2,1));
    G(1,:) = [];G(1,1) = -1/hs; G(end,end+1) = 1/hs;
    %%G(1,:) = [];G(1,1) = -1*0.5/hs; G(end,end+1) = 1*0.5/hs;
    A = -[sparse(m-2,m-2),G;-G',sparse(m,m)];
    A(1,m-1) = 0.5*A(1,m-1); A(m-2,end) = 0.5*A(m-2,end);
    
    %G = 1/(2*hs)*gallery('tridiag',-ones(m-2,1),zeros(m-1,1),ones(m-2,1));
    %G(1,:) = [];%G(1,1) = -1/hs;
    %G(end,end+1) = 0.5/hs;
    %A = -[sparse(m-1,m-1),G;-G',sparse(m-1,m-1)];
    %A(1,m-1) = 0.5*A(1,m-1); A(m-2,end) = 0.5*A(m-2,end);
    
    %cols = 1/(2*hs)*[-ones(m-2,1),sparse(m-2,1),ones(m-2,1)];
    %G = sparse(1:m-2,1:m-2,cols,m-2,m-1);
    %spdiags([sparse(m-2,1),sparse(m-2,1),-ones(m-2,1),sparse(m-2,1),ones()])
    
    
    
    %A = [sparse(m-2,m-2),G;G,sparse(m,m)];
    %G(1,2) = -1/hs; G(end,end-1) = -1/hs;
    %A = [sparse(m-2,m-2),G;G,sparse(m-2,m-2)];
    %J = [sparse(m,m-2),speye(m);-speye(m-2),sparse(m-2,m)];
    %J = [sparse(m-1,m-1),speye(m-1);-speye(m-1),sparse(m-1,m-1)];
    %J = [sparse(m,m),speye(m);-speye(m),sparse(m,m)];
    %isequal(J,-(J)')
    %isequal(J*A,(J*A)')
    %a = 2;
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
elseif strcmp(eqn,'maxwell3D') % DEnne MÅ testes!!!!!!!!!!
    %%%% ddx
    T = 1/(2*hs)*gallery('tridiag',-ones(m-1,1),zeros(m,1),ones(m-1,1));
    ddx = spalloc(m^3,m^3,4*m^2); % allokerer sikket for feil antall elementer
    for i = 1:m^2
        ddx(1+(i-1)*m:i*m,1+(i-1)*m:i*m) = T;
    end
    %%%% ddy
    ddy = spalloc(m^3,m^3,4*m^2);
    for i = 1:m^2
        ddy(m*(i-1)+1+1,2*m+i) = 1;
        ddy(m*(i-1)+1+1,i) = -1;
        ddy(m*(i-1)+1,m+i) = 1;
        ddy(m*i,m+i) = -1;
    end
    %%%% ddz
    ddz = spalloc(m^3,m^3,4*m^2);
    for i = 1:m^2
        ddz(m*(i-1)+1+1,2*m^2+i) = 1;
        ddz(m*(i-1)+1+1,i) = -1;
        ddz(m*(i-1)+1,m^2+i) = 1;
        ddz(m*i,m^2+i) = -1;
    end
    
    %%%% Making curl!
    my = 1; epsilon = 1; sigma = zeros(3*m^3,3*m^3);
    T = [sparse(m^3,m^3),-ddz,ddy;
        ddz,sparse(m^3,m^3),-ddx;
        -ddy,ddx,sparse(m^3,m^3)];
    
    A = [sparse(3*m^3,3*m^3), 1/my*T;1/epsilon*T,sigma];
elseif strcmp(eqn,'random')
    A = rand(2*(m-2)^2);
    A = 0.5*[sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)]*(A+A'+m^2*speye(2*(m-2)^2)); % eventuelt minus dette, men det har ikke noe å si!
elseif strcmp(eqn,'semirandom')
    try
        load('semirandomA.mat','A');
    catch
        A = -1;
    end
    if size(A,1) ~= 2*(m-2)^2
        D = abs(rand(2*(m-2)^2,1)) + 5*ones(2*(m-2)^2,1);
        D1 = rand(2*(m-2)^2-1,1);%rand(2*(m-2)^2-1,1);
        Atilde = gallery('tridiag',D1,D,D1);
        A = [sparse((m-2)^2,(m-2)^2),speye((m-2)^2);-speye((m-2)^2),sparse((m-2)^2,(m-2)^2)]*Atilde;
        save('semirandomA.mat','A');
    end
elseif strcmp(eqn,'eigen')
    
    try
        load('eigenQ.mat','Q');
    catch
        Q = -1;
    end
    
    if size(Q,1) ~= 2*(m-2)^2
        Q = rand(2*(m-2)^2);
        [Q,~] = qr(Q);
        save('eigenQ.mat','Q');
    end
    
    eigenvals = eigenvalues(2*(m-2)^2-1);
    D = gallery('tridiag',-eigenvals,zeros(2*(m-2)^2,1),eigenvals);
    %[a3,b3] = eigs(D);
    A =  Q*D*Q';
    %A1 =  Q*(expm(D))*Q';
    %[a,b] = eig(A,'vector');
    %[a1,b1] = eig(A1,'vector');
end
end

