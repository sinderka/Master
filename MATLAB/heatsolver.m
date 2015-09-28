function [ utdata ] = heatsolver(m,n,k,prob,solmeth,conv,para)


%%% Initsiell data
utdata = zeros(1,3);
Y = linspace(0,1,m+2)'; X = Y;
hs = 1/(m+1);
A = -1/hs^2*gallery('poisson', m);
T = linspace(0,1,k); ht = T(2)-T(1);
U = 0;
ant = 1;
%%%

    [f1,g1,f2,g2,solution] = getTestFunction('heat',prob);
    
if prob < 3 && solmeth < 3
    %[f1,g1,f2,g2,solution] = getTestFunction('heat',prob);
    v1 = zeros(m^2,1);
    v2 = zeros(m^2,1);
    for i = 1:m
        v1(1+(i-1)*m:i*m) = g1(X(i+1),Y(2:m+1));
        v2(1+(i-1)*m:i*m) = g2(X(i+1),Y(2:m+1));
    end
    F1 = zeros(1,k);
    F2 = zeros(1,k);
    for i = 1:k
        F1(i) = f1(T(i));
        F2(i) = f2(T(i));
    end
    if solmeth == 1 %% Restarted krylov
        tic;
        [U1,temp1] = KPMheat(F1,A,v1,k,ht,n,conv);
        [U2,temp2] = KPMheat(F2,A,v2,k,ht,n,conv);
                            %(Zn,A,v,k,ht,n,conv)
        U = U1+U2;
        utdata(1) = max(temp1,temp2);
        utdata(2) = toc;
    elseif solmeth == 2 %% Full krylov
        tic;
        [U1,~] = KPMheat(F1,A,v1,k,ht,m^2,1);
        [U2,~] = KPMheat(F2,A,v2,k,ht,m^2,1);
        U = U1+U2;
        utdata(1) = 1;
        utdata(2) = toc;
    end
else
    %delete(gcp)
    if prob < 3
        p1 =@(t,x,y) f1(t).*g1(x,y)+f2(t).*g2(x,y);
    else
        p1 = @(t,x,y) f1(t,x,y);
    end
    
    F1 = zeros(m^2,k);
    for j = 1:k
        for i = 1:m
            F1(1+(i-1)*m:i*m,j) = p1(T(j),X(i+1),Y(2:m+1));
        end
    end
    if solmeth == 1 %% Restarted krylov
        parpool(para);
        tic;
        parfor i = 1:m^2
            e_i = zeros(m^2,1); e_i(i) = 1;
            [U1,ant1] = KPMheat(F1(i,:),A,e_i,k,ht,n,conv);
            U = U + U1;
            ant = max(ant1,ant);
        end
        utdata(1) = ant;
        utdata(2) = toc;
        delete(gcp);
    elseif solmeth == 2 %% Full krylov
        tic;
        parpool(para);
        parfor i = 1:m^2
            e_i = zeros(m^2,1); e_i(i) = 1;
            U = U + KPMheat(F1(i,:),A,e_i,k,ht,m^2);
        end
        utdata(1) = 1;
        utdata(2) = toc;
        delete(gcp);
    elseif solmeth == 3 %% Direkte integrasjon
        tic;
        U = integrate(A,F1,m^2,k,ht,1);
        %(H,F,n,k,ht,hn)
        
        utdata(1) = 1;
        utdata(2) = toc;
    end
end

correctsolution = zeros(m^2,k);
for i = 1:k
    for j = 1:m
        correctsolution(1+(j-1)*(m):j*(m),i) = solution(T(i),X(j+1),Y(2:end-1));
    end
end

%%% Estimating error
Err = zeros(1,m^2);
for i = 1:m^2
    Err(i) = norm(correctsolution(i,:)-U(i,:) ,Inf)/max(abs(correctsolution(i,:)));
end
utdata(3) = max(Err);
%%%