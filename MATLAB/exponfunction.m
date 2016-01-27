function exponmat = exponfunction(A)
exponmat = eye(size(A));
for i = 1:100
    exponmat = exponmat + A^i/factorial(i);
end


end

