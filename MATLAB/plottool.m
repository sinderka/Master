function plottool(m,n,k,eqn,alg,restart,prob,conv,para,leg,lab,data,type,help,name)
%m,n,k,para,conv kan være lister med tall
%eqn, alg, restart,prob, kan ikke være lister med tall
%help sier noe om hvordan hjelpelinjen skal se ut
%type sier noe om hvilke type plott som skal benyttes
%data sier noe om hva av utdataen som skal plottes
%leg, lab og name angir tekst og navn på figuren

linetype = {'k:+','k:o','k:*','k:.','k:x','k:s','k:d','k:^','k:v','k:<','k:>','k:p','k:h'};

%%%%% TODO %%%%
% alg og restart burde kunne være lister!
% Automatisk legge til akser og navn for å lettere forhindre feil?


[p,b] = getPandB(m,n,k,para,conv);
ant1 = length(p)-1; ant2 = length(b)-1;
a = 1; b = 1; c = 1; d = 1; e = 1;
utdata = zeros(ant2,ant1,5);
for i = 1:ant1
    [a0,b0,c0,d0,e0] = addOne(m,n,k,para,conv,a,c,b,d,e,-1);
    a = max(1,a0*(a0-a)); b = max(1,b0*(b0-b)); c = max(1,c0*(c0-c)); d = max(1,d0*(d0-d)); e = max(1,e0*(e0-e));
    for j = 1:ant2
        [a,b,c,d,e] = addOne(m,n,k,para,conv,a,c,b,d,e,-2);
        utdata(j,i,:) = solver(m(a),n(b),k(c),eqn,alg,restart,prob,conv(d),para(e));
    end
end

for i = 1:ant2
    if strcmp(type,'plot')
        plot(p(2:end),utdata(:,data,i),char(linetype(i)))
    elseif strcmp(type,'loglog')
        loglog(p(2:end),utdata(:,data,i),char(linetype(i)))
    elseif strcmp(type,'semilogx')
        semilogx(p(2:end),utdata(:,data,i),char(linetype(i)))
    elseif strcmp(type,'semilogy')
        semilogy(p(2:end),utdata(i,:,data),char(linetype(i)))
    end
    hold on
end

if help(1)
    plot(p,p.^help(2),'k-')
    leg = {leg , 'Helpline'};
end



ylabel(lab(1));
xlabel(lab(2));
legend(char(leg));
h = set(findall(gcf,'-property','FontSize'), 'Fontsize',18);
set(h,'Location','Best');
location = strcat('/home/shomeb/s/sindreka/Master/MATLAB/fig/',name);
saveas(gcf,location,'fig');
saveas(gcf,location,'jpeg');
hold off
end

function [p,b] = getPandB(m,n,k,para,conv)
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
    
    if length(conv) > 1
        if conv(1) == -1
            p = conv;
        elseif conv(1) == -2
            b = conv;
        end
    end
end
function [a,b,c,d,e] = addOne(m,n,k,para,conv,a,c,b,d,e,var)
if m(1) == var
    a = a + 1;
elseif n(1) == var
    b = b + 1;
elseif k(1) == var
    c = c + 1;
elseif para(1) == var
    d = d + 1;
elseif conv(1) == var
    e = e + 1;
end
end