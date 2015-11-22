function savedatastorage(  m,n,k,restart,eqn,prob,para,alg,conv ,data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if alg < 3
    str1 = {num2str(m),num2str(n),num2str(k),num2str(restart),num2str(eqn),num2str(prob),num2str(para),num2str(alg),num2str(conv)};
    try
        datalist = load('saveddata12.mat');
    catch
        temp = [str1,data];
        save('saveddata12.mat','temp')
        return
    end
    datalist(end+1,:) = [str1,data];
    save('saveddata12.mat',datalist)
    return
end
str1 = {num2str(m),num2str(k),num2str(eqn),num2str(prob),num2str(para)};
try
    datalist = load('saveddata3.mat');
catch
    temp = [str1,data];
    save('saveddata3.mat','temp')
    return
end
datalist(end+1,:) = [str1,data];
save('saveddata3.mat',datalist)
return
end

