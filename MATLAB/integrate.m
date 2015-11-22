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
%mat = inv(speye(n)-A*ht/2); U(:,2) = mat*ht/2*(F(:,2)+F(:,1));


U(:,2) = ht/2 * (F(:,1) + F(:,2) + ht*A*F(:,1) );

for j = 3:1:k
    %U(:,j) = mat*(U(:,j-1) + ht/2*A*U(:,j-1)+ht/2*(F(:,j)+F(:,j-1)));
    
    U(:,j) = U(:,j-1) +ht/2*A*U(:,j-1) + ht/2*F(:,j-1) + ht/2*A*U(:,j-1)+ht^2/2*A^2*U(:,j-1)  + ht^2/2*A*F(:,j-1) + ht/2*F(:,j);
    
end

end
