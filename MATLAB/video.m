function video(U,m,T,eqn)
figure(1)
if max(max(U)) == min(min(U))
    display('Something went wrong with the video!')
    return
end
k = size(U,2);

if strcmp(eqn,'maxwell1D')
    for i = 1:k
        plot(U(:,i))
        caxis([min(min(U)),max(max(U))])
        axis([0,m,min(min(U)),max(max(U))])
        drawnow
        pause(T)
    end
    return
end

for i = 1:k
    mesh(reshape(U(:,i),m,m))
    caxis([min(min(U)),max(max(U))])
    axis([0,m,0,m,min(min(U)),max(max(U))])
    drawnow
    pause(T)
end

end

