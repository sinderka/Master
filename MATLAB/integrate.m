% Må nok testes litt
function Zn = integrate(H,F,n,k,ht) %%% Intigrere direkte
%%%%%Solves the equation z'-H*z= F Numerically
%%% ht is stepsize, n is the size of the square matrix H
%%% k is the number of time-steps.
%%% Zn(:,1) = 0.
if min(size(F)) == 1
    
    mat = inv(eye(n)-ht/2*H);
    Zn = zeros(n,k);
    Zn(:,2) = mat*(ht/2*(F(:)+F(:)));
    for j = 3:1:k
        Zn(:,j) = mat*(Zn(:,j-1)+ht/2*(H*Zn(:,j-1)+(F(:)+F(:))));
    end
else
    %function Zn = locintegrate2(H,F,n,k,ht)
    %%%%%Solves the equation z'-H*z=F Numerically
    %%% ht is stepsize, n is the size of the square matrix H
    %%% k is the number of time-steps.
    
    mat = inv(eye(n)-ht/2*H);
    Zn = zeros(n,k);
    Zn(:,2) = mat*(ht/2*(F(:,1)+F(:,2)));
    for j = 3:1:k
        Zn(:,j) = mat*(Zn(:,j-1)+ht/2*(H*Zn(:,j-1)+(F(:,j)+F(:,j-1))));
    end
end
end