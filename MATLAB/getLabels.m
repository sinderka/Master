function [ylab,xlab,leg,additionalInfo] = getLabels(ant2,m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,PMint)
if data == 1
    ylab = {'r_n'};
elseif data == 2
    ylab = {'T_c'};
elseif data == 3
    ylab = {'er_1'};
elseif data == 4
    ylab = {'en_1'};
elseif data == 5
    ylab = {'er_2'};
elseif data == 6
    ylab = {'en_2'};
elseif data == 7
    ylab = {'en_3'};
elseif data == 8
    ylab = {'en_4'};
end
if m(1) == -1
    xlab = {'m'};
end
if n(1) == -1
    xlab = {'n'};
end
if simtime(1) == -1
    xlab = {'T_s'};
end
if K(1) == -1
    xlab = {'K'};
end
if k(1) == -1
    xlab = {'k'};
end
if para(1) == -1
    xlab = {'p_n'};
end
if conv(1) == -1
    xlab = {'\epsilon'};
end
if alg(1) == -1
    xlab = {'solution method'};
end
if int(1) == -1
    xlab = {'integration method'};
end
if restart(1) == -1
    xlab = {'r_n'};
end
if prob(1) == -1
    xlab = {'problem'};
end
if PMint(1) == -1
    xlab{'PM integration method'};
end
if  ~exist('xlab','var')
    xlab = {'T_s'};
end
if m(1) == -1 && k(1) == -1
    if isequal(m,k)
        xlab = {'m=k'};
    else
        xlab = {'m^2 = k'};
    end
end
if K(1) == -1 && k(1) == -1
    stri = num2str(max(K.*k));
    stri = strcat('K \cdot k = ',stri);
    xlab = {stri};
end
if k(1) == -1 && simtime(1) == -1
    stri = num2str(max(k./simtime));
    stri = strcat(stri,'\cdot T_s', '=k');
    xlab = {stri};
end

leg = {};
if m(1) == -2
    for i = 1:ant2
        stri = strcat('m=',num2str(m(i+1)));
        leg(i) = {stri};
    end
elseif n(1) == -2
    for i = 1:ant2
        stri = strcat('n=',num2str(n(i+1)));
        leg(i) = {stri};
    end
elseif simtime(1) == -2
    for i = 1:ant2
        stri = strcat('T_s: ',num2str(simtime(i+1)));
        leg(i) = {stri};
    end
elseif K(1) == -2
    for i = 1:ant2
        stri = strcat('K: ',num2str(K(i+1)));
        leg(i) = {stri};
    end
elseif k(1) == -2
    for i = 1:ant2
        stri = strcat('k=',num2str(k(i+1)));
        leg(i) = {stri};
    end
elseif para(1) == -2
    for i = 1:ant2
        stri = strcat('p_n=',num2str(para(i+1)));
        leg(i) = {stri};
    end
elseif conv(1) == -2
    for i = 1:ant2
        stri = strcat('\epsilon=',num2str(conv(i+1)));
        leg(i) = {stri};
    end
elseif alg(1) == -2
    for i = 1:ant2
        if alg(i+1) == 1
            stri = 'KPM';
        elseif alg(i+1) == 2
            stri = 'SLPM';
        elseif alg(i+1) == 3
            stri = 'DM';
        end
        leg(i) = {stri};
    end
elseif int(1) == -2
    for i = 1:ant2
        if int(i+1) == 1
            stri = 'trapezoidal rule';
        elseif int(i+1) == 2
            stri = 'forward Euler';
        elseif int(i+1) == 3
            stri = 'midpoint rule';
        end
        leg(i) = {stri};
    end
elseif restart(1) == -2
    for i = 1:ant2
        stri = strcat('r_n=',num2str(restart(i+1)));
        leg(i) = {stri};
    end
elseif prob(1) == -2
    for i = 1:ant2
        stri = strcat('problem=',num2str(prob(i+1)));
        leg(i) = {stri};
    end
elseif PMint(1) == -2
    for i = 1:ant2
        if PMint(i+1) == 1
            if int == 1
                stri = 'trapezoidal rule';
            elseif int == 2
                stri = 'forward Euler';
            elseif int == 3
                stri = 'midpoint rule';
            end
        elseif PMint(i+1) == 2
            stri = 'Matlab expm';
        elseif PMint(i+1) == 3
            stri = 'My expm';
        end
        leg(i) = {stri};
    end
end
additionalInfo = {eqn};

if length(m) == 1
    stri = strcat('m=',num2str(m));
    additionalInfo(end+1) = {stri};
end
if length(n) == 1
    stri = strcat('n=',num2str(n));
    additionalInfo(end+1) = {stri};
end
if length(simtime) == 1
    stri = strcat('T_s: ',num2str(simtime));
    additionalInfo(end+1) = {stri};
end
if length(K) == 1
    stri = strcat('K=',num2str(K));
    additionalInfo(end+1) = {stri};
end
if length(k) == 1
    stri = strcat('k=',num2str(k));
    additionalInfo(end+1) = {stri};
end
if length(alg) == 1
    if alg == 1
        stri = 'KPM';
    elseif alg == 2
        stri = 'SLPM';
    elseif alg == 3
        stri = 'DM';
    end
    additionalInfo(end+1) = {stri};
end

if length(int) == 1
    if int == 1
        stri = 'trapezoidal rule';
    elseif int == 2
        stri = 'forward Euler';
    elseif int == 3
        stri = 'midpoint rule';
    end
    additionalInfo(end+1) = {stri};
end
if length(restart) == 1
    stri = strcat('restart=',num2str(restart));
    additionalInfo(end+1) = {stri};
end
if length(prob) == 1
    stri = strcat('problem=',num2str(prob));
    additionalInfo(end+1) = {stri};
end
if length(conv) == 1
    stri = strcat('\epsilon=1e',num2str(log10(conv)));
    additionalInfo(end+1) = {stri};
end
if length(para) == 1
    stri = strcat('p_n=',num2str(para));
    additionalInfo(end+1) = {stri};
end
if length(PMint) == 1
    if PMint == 1
        if int == 1
            stri = 'trapezoidal rule';
        elseif int == 2
            stri = 'forward Euler';
        elseif int == 3
            stri = 'midpoint rule';
        end
    elseif PMint == 2
        stri = 'Matlabs built in expm';
    elseif PMint == 3
        stri = 'eigenvalue + exp';
    end
    additionalInfo(end+1) = {stri};
end
end