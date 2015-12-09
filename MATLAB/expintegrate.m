function U = expintegrate(A,F,k,ht)
%solves the problem du/dt = Au+F
%indata
% A: mxm matrix
% F: k row
% k: number of points in time
% ht: step size in time
%outdata
% U: the solution

U = zeros(size(F));
for i = 0:k-1
    U(:,i+1) = expm(A*i*ht)*F(:,i+1);
end
end
