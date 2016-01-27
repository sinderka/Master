function eigenvals = eigenvalues(m)

eigenvals = zeros(1,m);

for j = 2:m+1
    eigenvals(j-1) = -2*j^2/2^j;
end

end