% Må nok testes litt
function Zn = test2(H,F,n,k,ht,hn) %%% Intigrere direkte
%%%%%Solves the equation z'-H*z=F Numerically
%%% ht is stepsize, n is the size of the square matrix H
%%% k is the number of time-steps.
%%% Zn(:,1) = 0.
    mat = inv(eye(n)-ht/2*H);
    %e1 = zeros(n,1); e1(1) = 1;
    e1 = 1;
    Zn = zeros(n,k);
    Zn(:,2) = mat*(ht/2*hn*e1*(F(:)+F(:)));
    for j = 3:1:k
        Zn(:,j) = mat*(Zn(:,j-1)+ht/2*(H*Zn(:,j-1)+hn*e1*(F(:)+F(:))));
    end
end