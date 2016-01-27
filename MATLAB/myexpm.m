function U = myexpm(A,b,T)

U = zeros(length(A),length(T));
[V,D] = eig(A);
for i = 1:length(T)
U(:,i) = V*diag(exp(diag(D*T(i))))/V * b - b;
end

end