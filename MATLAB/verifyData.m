function verified = verifyData(m,n,k,eqn,alg,restart,prob,para)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%%%% Denne må stadig oppdateres!

if (m > 0 && k > 0) && (n > 0 && prob > 0) && alg > 0 && (restart == 0 || restart == 1)
    
    if (strcmp(eqn,'heat') && prob < 6 && para <= 8 && (alg == 1 || alg == 3 ) && n <= (m-2)^2)
        verified = 1;
        return
    elseif (strcmp(eqn,'wave') && prob < 7 && para <= 8 && (alg == 1 || alg == 2 || alg == 3) && n <= (m-2)^2)
        verified = 1;
        return
    elseif (strcmp(eqn,'maxwell') && prob < 0 && para <= 8 && (alg == 1 || alg == 2 || alg == 3) && n <= (m-2)^3) % Fyll inn!
        verified = 1;
        return
    end
    
end
verified = 0;
end

