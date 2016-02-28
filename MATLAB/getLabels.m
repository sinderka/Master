function [ylab,xlab,leg,additionalInfo,helpinfo] = getLabels(ant2,m,n,simtime,K,k,eqn,alg,int,restart,prob,conv,para,data,PMint,iter)
helpinfo = -1;
if data(1) == 1
    ylab = {'\it i_r'};
elseif data(1) == 2
    ylab = {'\it T_c'};
elseif data(1) == 3
    ylab = {'\it error'};
elseif data(1) == 4
    ylab = {'energy \it H_1'};
elseif data(1) == 5
    ylab = {'\it error^{\rm comp}'};
elseif data(1) == 6
    ylab = {'energy \it H_1^{\rm comp}'};
elseif data(1) == 7
    ylab = {'energy \it H_3'};
elseif data(1) == 8
    ylab = {'energy \it H_4'};
elseif data(1) == 9
    ylab = {'abs(\it H_3-H_4\rm)'};
else
    ylab = {'energy \it H_1'};
end
if m(1) == -1
    xlab = {'\it m'};
    helpinfo = '\it m';
end
if n(1) == -1
    xlab = {'\it n'};
    helpinfo = '\it n';
end
if simtime(1) == -1
    xlab = {'\it T_s'};
end
if K(1) == -1
    xlab = {'\it K'};
    helpinfo = '\it K';
end
if k(1) == -1
    xlab = {'\it k'};
    helpinfo = '\it k';
end
if para(1) == -1
    xlab = {'p_n'};
end
if conv(1) == -1
    if ~(min(conv(2:end)) > 0.99)
        xlab = {'\iota'};
    else
        xlab = {'i_r'};
    end
end
if alg(1) == -1
    xlab = {'solution method'};
end
if int(1) == -1
    xlab = {'integration method'};
end
if restart(1) == -1
    xlab = {'\it i_n'};
end
if prob(1) == -1
    xlab = {'problem'};
end
if PMint(1) == -1
    xlab{'PM integration method'};
end
if  ~exist('xlab','var')
    xlab = {'\it T_s'};
end
if m(1) == -1 && k(1) == -1
    if isequal(m,k)
        xlab = {'\it k, with m = k'};
    else
        xlab = {'\it k\rm , with \it m^2 = k'};
    end
end
if K(1) == -1 && k(1) == -1
    stri = num2str(max(K.*k));
    stri = ['\it k \rm, with \it K \cdot k = ',stri];
    xlab = {stri};
end
if k(1) == -1 && simtime(1) == -1
    stri = num2str(max(k./simtime));
    stri = ['\it T_s \rm, with \it k = ',stri,' \cdot T_s'];
    xlab = {stri};
end
if K(1) == -1 && simtime(1) == -1
    if (max(K./simtime)) == 1
        stri = '\it T_s \rm, with \it K = T_s';
    else
        stri = num2str(max(K./simtime));
        stri = ['\it T_s \rm, with \it K = ',stri,' \cdot T_s'];
    end
    
    xlab = {stri};
end


leg = {};
if m(1) == -2
    for i = 1:ant2
        if length(n) == 1
            nstr = num2str(n);
        elseif isequal(m,n)
            nstr = 'm';
        else
            nstr = 'n';
        end
        if alg == 1
            tstr = ['KPM(',nstr,')'];
        elseif alg == 2
            tstr = ['SLM(',nstr,')'];
        elseif alg == 3 && ~( data == 1 || data == 5 || data == 6 )
            tstr = ['DM'];
        end
        
        stri = ['m=',num2str(m(i+1))];
        leg(i) = {stri};
    end
elseif data(1) == -3
    if length(n) == 1
        nstr = num2str(n);
    else
        nstr = 'n';
    end
    for i = 2:length(data)
        if data(i) == 13
            stri = 'DM';
        elseif data(i) == 10
            stri = ['z_{',nstr,'}(t)'];
        elseif data(i) == 11
            stri = ['u_{',nstr,'}(t)'];
            %elseif data(i) == 12
            %    stri = ['z_',nstr,'^{(' ,iter, ')}(t)'];
        elseif data(i) == 4
            stri = ['u_{',nstr,'}^{(' ,iter, ')}(t)'];
        elseif data(i) == 7
            stri = 'energy \it H_3';
        elseif data(i) == 8
            stri = 'energy \it H_4';
        elseif data(i) == 9
            stri = 'abs(\it H_3-H_4\rm)';
            ylab = {'energy'};
        end
        leg(i-1) = {stri};
    end
elseif n(1) == -2
    if length(alg) == 1
        if alg == 1
            tstr = 'KPM(';
        elseif alg == 2
            tstr = 'SLM(';
        end
        for i = 1:ant2
            stri = [tstr,num2str(n(i+1)),')'];
            leg(i) = {stri};
        end
    else
        for i = 1:ant2
            if alg(i+1) == 1
                tstr = 'KPM(';
            elseif alg(i+1) == 2
                tstr = 'SLM(';
            end
            stri = [tstr,num2str(n(i+1)),')'];
            leg(i) = {stri};
        end
    end
elseif simtime(1) == -2
    for i = 1:ant2
        stri = ['T_s = ',num2str(simtime(i+1))];
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
        stri = strcat('\iota=',num2str(conv(i+1)));
        leg(i) = {stri};
    end
elseif alg(1) == -2
    for i = 1:ant2
        if length(n) == 1;
            tstr = num2str(n);
        else
            tstr = 'n';
        end
        
        if alg(i+1) == 1
            if PMint(1) == -2
                if PMint(i+1) == 1
                    str = ' trap';
                elseif PMint(i+1) == 2
                    str = ' expm';
                elseif PMint(i+1) == 3
                    str = ' diag';
                end
            else
                str = '';
            end
            stri = ['KPM(',tstr,')' , str];
            leg(i) = {stri};
        elseif alg(i+1) == 2
            if PMint(1) == -2
                if PMint(i+1) == 1
                    str = ' trap';
                elseif PMint(i+1) == 2
                    str = ' expm';
                elseif PMint(i+1) == 3
                    str = ' diag';
                end
            else
                str = '';
            end
            stri = ['SLM(',tstr,')',str];
            leg(i) = {stri};
        elseif alg(i+1) == 3 && ~( data == 1 || data == 5 || data == 6 )
            stri = 'DM';
            leg(i) = {stri};
        end
        
    end
elseif int(1) == -2
    for i = 1:ant2
        if int(i+1) == 1
            stri = 'trap';
        elseif int(i+1) == 2
            stri = 'Euler';
        elseif int(i+1) == 3
            stri = 'mid';
        end
        leg(i) = {stri};
    end
elseif restart(1) == -2
    for i = 1:ant2
        stri = strcat('\it i_n=',num2str(restart(i+1)));
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
                stri = 'trap';
            elseif int == 2
                stri = 'Euler';
            elseif int == 3
                stri = 'mid';
            end
        elseif PMint(i+1) == 2
            stri = 'expm';
        elseif PMint(i+1) == 3
            stri = ' diag';
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
if length(alg) == 1 || length(alg) == 2
    if length(n) == 1;
        tstr = num2str(n);
    else
        tstr = 'n';
    end
    if alg(end) == 1
        stri = strcat('KPM(',tstr,')');
    elseif alg(end) == 2
        stri = strcat('SLM(',tstr,')');
    elseif alg(end) == 3 
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
    stri = strcat('\iota=1e',num2str(log10(conv)));
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
        stri = 'Matlabs expm function';
    elseif PMint == 3
        stri = 'Eigenvalue and diagonalization';
    end
    additionalInfo(end+1) = {stri};
end
end