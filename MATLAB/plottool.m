function plottool(m,n,k,eqn,alg,int,restart,prob,conv,para,data,type,help,name,save)
%m,n,k,para,conv,alg,restart,prob kan være lister med tall
%eqn Angir hvilken ligning som løses
%help sier noe om hvordan hjelpelinjen skal se ut
%type sier noe om hvilke type plott som skal benyttes
%data sier noe om hva av utdataen som skal plottes
%leg, lab og name angir tekst og navn på figuren

linetype = {'k:+','k:o','k:*','k:.','k:x','k:s','k:d','k:^','k:v','k:<','k:>','k:p','k:h'};




[p,b] = getPandB(m,n,k,alg,int,restart,prob,conv,para);
ant1 = length(p)-1; ant2 = length(b)-1;
a = 1; b = 1; c = 1; d = 1; e = 1; f = 1; g = 1; h = 1; aa = 1;
utdata = zeros(ant2,ant1,6);
for i = 1:ant1
    [a0,b0,c0,d0,e0,f0,g0,h0,aa0] = addOne(m,n,k,alg,int,restart,prob,conv,para,a,c,b,d,e,f,g,h,aa,-1);
    a = max(1,a0*(a0-a)); b = max(1,b0*(b0-b)); c = max(1,c0*(c0-c)); d = max(1,d0*(d0-d)); e = max(1,e0*(e0-e)); f = max(1,f0*(f0-f)); g = max(1,g0*(g0-g)); h = max(1,h0*(h0-h)); aa = max(1,aa0*(aa0-aa));
    for j = 1:ant2
        [a,b,c,d,e,f,g,h,aa] = addOne(m,n,k,alg,int,restart,prob,conv,para,a,c,b,d,e,f,g,h,aa,-2);
        utdata(j,i,:) = solver(m(a),n(b),k(c),eqn,alg(f),int(aa),restart(g),prob(h),conv(d),para(e));
        %utdata(j,i,:) = energyTest(m(a),n(b),k(c),eqn,alg(f),int(aa),restart(g),prob(h),conv(d),para(e));
    end
end

for i = 1:ant2
    if strcmp(type,'plot')
        plot(p(2:end),utdata(i,:,data),char(linetype(i)))
    elseif strcmp(type,'loglog')
        loglog(p(2:end),utdata(i,:,data),char(linetype(i)))
    elseif strcmp(type,'semilogx')
        semilogx(p(2:end),utdata(i,:,data),char(linetype(i)))
    elseif strcmp(type,'semilogy')
        semilogy(p(2:end),utdata(i,:,data),char(linetype(i)))
    end
    hold on
end


[ylab,xlab,leg,additionalInfo] = getLabels(ant2,m,n,k,eqn,alg,int,restart,prob,conv,para,data);
%Må lage kode for legende!!!! Må jeg lage noe for tittel?
%text(0.25,2.5,xlab,'Interpreter','latex')
if help(1)
    plot(p(2:end),p(2:end).^help(2),'k-')
    leg(end+1) = {'Helpline'};
end
ylabel(ylab);
xlabel(xlab);
%set(gcf,'CurrentAxes',h)
%text(1,1,char(additionalInfo),'FontSize',12)
title(char(additionalInfo))
%dim = [0.9 0.1 0.3 0.3];
%annotation('textbox',dim,'String',char(additionalInfo),'FitBoxToText','on');

legend(char(leg));
h = set(findall(gcf,'-property','FontSize'), 'Fontsize',18);
set(h,'Location','Best');
if save
    location = strcat('/home/shomeb/s/sindreka/Master/MATLAB/fig/',name);
    saveas(gcf,location,'fig');
    saveas(gcf,location,'jpeg');
end
hold off
end

function [p,b] = getPandB(m,n,k,alg,int,restart,prob,conv,para)
if length(m) > 1
    if m(1) == -1
        p = m;
    elseif m(1) == -2
        b = m;
    end
end

if length(n) > 1
    if n(1) == -1
        p = n;
    elseif n(1) == -2
        b = n;
    end
end

if length(k) > 1
    if k(1) == -1
        p = k;
    elseif k(1) == -2
        b = k;
    end
end

if length(para) > 1
    if para(1) == -1
        p = para;
    elseif para(1) == -2
        b = para;
    end
end
if length(alg) > 1
    if alg(1) == -1
        p = alg;
    elseif alg(1) == -2
        b = alg;
    end
end
if length(int) > 1
    if int(1) == -1
        p = int;
    elseif int(1) == -2
        b = int;
    end
end
if length(restart) > 1
    if restart(1) == -1
        p = restart;
    elseif restart(1) == -2
        b = restart;
    end
end
if length(prob) > 1
    if prob(1) == -1
        p = prob;
    elseif prob(1) == -2
        b = prob;
    end
end
if length(conv) > 1
    if conv(1) == -1
        p = conv;
    elseif conv(1) == -2
        b = conv;
    end
end
end

function [a,b,c,d,e,f,g,h,aa] = addOne(m,n,k,alg,int,restart,prob,conv,para,a,c,b,d,e,f,g,h,aa,var)
if m(1) == var
    a = a + 1;
end
if n(1) == var
    b = b + 1;
end
if k(1) == var
    c = c + 1;
end
if para(1) == var
    d = d + 1;
end
if conv(1) == var
    e = e + 1;
end
if alg(1) == var
    f = f + 1;
end
if restart(1) == var
    g = g + 1;
end
if prob(1) == var
    h = h + 1;
end
if int(1) == var
    aa = aa + 1;
end
end

function [ylab,xlab,leg,additionalInfo] = getLabels(ant2,m,n,k,eqn,alg,int,restart,prob,conv,para,data)
if data == 1
    ylab = {'iterations'};
elseif data == 2
    ylab = {'time'};
elseif data == 3
    ylab = {'error'};
elseif data == 4
    ylab = {'energy error'};
elseif data == 5
    ylab = {'difference error'};
elseif data == 6
    ylab = {'difference energy'};
end
if m(1) == -1
    xlab = {'m'};
elseif n(1) == -1
    xlab = {'n'};
elseif k(1) == -1
    xlab = {'k'};
elseif para(1) == -1
    xlab = {'# processors'};
elseif conv(1) == -1
    xlab = {'convergence criterion'};
elseif alg(1) == -1
    xlab = {'solution method'};
elseif int(1) == -1
    xlab = {'solution method'};
elseif restart(1) == -1
    xlab = {'restart'};
elseif prob(1) == -1
    xlab = {'problem'};
end
leg = {};
if m(1) == -2
    for i = 1:ant2
        stri = strcat('m=',num2str(m(i+1)));
        leg(i) = {stri};
    end
elseif n(1) == -2
    for i = 1:ant2
        stri = strcat('n=',num2str(m(i+1)));
        leg(i) = {stri};
    end
elseif k(1) == -2
    for i = 1:ant2
        stri = strcat('k=',num2str(m(i+1)));
        leg(i) = {stri};
    end
elseif para(1) == -2
    for i = 1:ant2
        stri = strcat('# processors=',num2str(m(i+1)));
        leg(i) = {stri};
    end
elseif conv(1) == -2
    for i = 1:ant2
        stri = strcat('convergence criterion=',num2str(m(i+1)));
        leg(i) = {stri};
    end
elseif alg(1) == -2
    for i = 1:ant2
        if alg(i+1) == 1
            stri = 'Arnoldi';
        elseif alg(i+1) == 2
            stri = 'SympLan';
        elseif alg(i+1) == 3
            stri = 'DirMet';
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
        stri = strcat('restart=',num2str(restart(i+1)));
        leg(i) = {stri};
    end
elseif prob(1) == -2
    for i = 1:ant2
        stri = strcat('problem=',num2str(prob(i+1)));
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
if length(k) == 1
    stri = strcat('k=',num2str(k));
    additionalInfo(end+1) = {stri};
end
if length(alg) == 1
    if alg == 1
        stri = 'Arnoldi';
    elseif alg == 2
        stri = 'SympLan';
    elseif alg == 3
        stri = 'DirMet';
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
    stri = strcat('prob=',num2str(prob));
    additionalInfo(end+1) = {stri};
end
if length(conv) == 1
    stri = strcat('convergence criterion=1e',num2str(log10(conv)));
    additionalInfo(end+1) = {stri};
end
if length(para) == 1
    stri = strcat('# processors=',num2str(para));
    additionalInfo(end+1) = {stri};
end
end