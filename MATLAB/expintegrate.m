function U = expintegrate(A,F,k,ht)
if size(F,2) == 1
    U = zeros(length(A),k);
    for i = 0:k-1
        U(:,i+1) = expm(A*i*ht)*F;
    end
else
    U = zeros(length(F),length(T));
    for i = 0:k-1
        U(:,i+1) = expm(A*i*ht)*F(:,i+1);
    end
end
end
