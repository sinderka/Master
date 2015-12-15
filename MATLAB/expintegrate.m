function U = expintegrate(A,F,k,ht)
%solves the problem du/dt = Au+F
%indata
% A: mxm matrix
% F: k row
% k: number of points in time
% ht: step size in time
%outdata
% U: the solution

if size(F,2) == 1
    %U = zeros(size(F));
    U = zeros(length(A),k);
    for i = 0:k-1
    %U(:,i+1) = expm(A*i*ht)*F(:,i+1);
        U(:,i+1) = expm(A*i*ht)*F;
    end
else
    %U = zeros(size(F,1),size(F,2));
    U = zeros(size(F));
    for i = 0:k-1
        U(:,i+1) = expm(A*i*ht)*F(:,i+1);
    end
 end
end