function [ bool1,data ] = datastorage( m,n,k,restart,eqn,prob,para,alg,conv )
bool1 = 0;
data = zeros(1,5);
if alg < 3
    str1 = {num2str(m),num2str(n),num2str(k),num2str(restart),num2str(eqn),num2str(prob),num2str(para),num2str(alg),num2str(conv)};
    try
        datalist = load('saveddata12.mat');
    catch
        bool1 = 0;
        data = zeros(1,5);
        return
    end
    for i = 1:size(datalist,1) % Er det problematisk at str listen kan ha forskjellige lengder?
        if strcmp(str1, datalist(i,1))
            bool1 = 1;
            data = datalist(i,2:end);
        end
    end
    return
end
str1 = {num2str(m),num2str(k),num2str(eqn),num2str(prob),num2str(para)};
try
    datalist = load('saveddata3.mat');
catch
    bool1 = 0;
    data = zeros(1,5);
    return
end
for i = 1:size(datalist,1) % Er det problematisk at str listen kan ha forskjellige lengder?
    if strcmp(str1, datalist(i,1))
        bool1 = 1;
        data = datalist(i,2:end);
    end
end

