% function U = integrate(A,F,k,ht)
% %Solves the equation u'-A*u= F numerically with the trapezoidal rule
% % A : A square matrix
% % ht: Stepsize
% % F : Left hand function
% %Returns:
% % U : A solution to the problem
% 
% n = size(A,1);
% U = zeros(n,k);
% mat = inv(speye(n)-(ht/2)*A);
% 
% %if size(F,2) == 1
% %    U(:,2) = mat*(ht*F);
% %    for j = 3:1:k
% %        U(:,j) = mat*(U(:,j-1)+ht/2*(A*U(:,j-1)+2*F));
% %    end
% %else
%     U(:,2) = mat*(ht/2*(F(:,1)+F(:,2)));
%     for j = 3:1:k
%         U(:,j) = mat*(U(:,j-1)+ht/2*(A*U(:,j-1)+(F(:,j)+F(:,j-1))));
%     end
% %end
% end

function U = integrate(A,F,k,ht)
%Solves the equation u'-A*u= F numerically with the midpoint rule
n = size(A,1);
U = zeros(n,k);

% Trapez method
mat = inv(speye(n)-A*ht/2); %U(:,2) = mat*ht/2*(F(:,2)+F(:,1));

% ""Modified Euler""
%U(:,2) = ht/2 * (F(:,1) + F(:,2) + ht*A*F(:,1) );

% ""Forward Euler""
%U(:,2) = ht*(A*U(:,1)+F(:,1));
%ht = 1/k;
% Implisitt Midpoint rule dobbel steglengde
%U(:,2) = ( speye(n) - A*ht )\(2*ht * (F(:,1) ));
% Eksplisitt Midpoint rule dobbel steglengde
%U(:,2) = 2*ht*( F(:,1) ); %U(:,3) = 2*ht*F(:,2); %
%U(:,2) = U(:,1) + ht*( A*( U(:,j-1) + ht/2*( A*U(:,j-1) + F(:,j-1) ) ) + 1/2*( F(:,j) + F(:,j-1) ) );
% Midpoint rule approximert F
%U(:,2) =
%ht = 1/k;
for j = 2:1:k
% Trapez method
    U(:,j) = mat*(U(:,j-1) + ht/2*A*U(:,j-1)+ht/2*(F(:,j)+F(:,j-1)));
    
% ""Modified Euler""
%    U(:,j) = U(:,j-1) +ht/2*A*U(:,j-1) + ht/2*F(:,j-1) + ht/2*A*U(:,j-1)+ht^2/2*A^2*U(:,j-1)  + ht^2/2*A*F(:,j-1) + ht/2*F(:,j);
    
% ""Forward Euler""
%    U(:,j) = U(:,j-1) + ht*(A*U(:,j-1)+F(:,j-1));
        
% Eksplisitt midpoint rule dobbel steglengde % Energien er periodisk
%    U(:,j) = U(:,j-2) + 2*ht * ( A*( U(:,j-2)+ht*( A*U(:,j-2)+F(:,j-2) ) ) + F(:,j-1) );

% Eksplisitt midpoint rule approximert F % Fungerer ikke !! % Energien
% stiger lineært
%   U(:,j) = U(:,j-1) + ht*( A*( U(:,j-1) + ht/2*( A*U(:,j-1) + F(:,j-1) ) ) + F(:,j) );

% Implisitt midpoint rule dobbel steglengde
%    U(:,j) = ( speye(n) - A*ht )\(U(:,j-2) + 2*ht * ( A*U(:,j-2)/2 + F(:,j-1) ));

% Implisitt midpoint rule approximert F % Fungerer!!
%    U(:,j) = mat*(U(:,j-1) + ht/2*( A*U(:,j-1) + F(:,j) + F(:,j-1) ));

end

end
