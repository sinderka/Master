function video(U,m,k,T)

%%% TODO %%%%
% Fikse fargene på videoen

for i = 1:k
    mesh(reshape(U(:,i),m,m))
    caxis([min(min(U)),max(max(U))])
    axis([0,m,0,m,min(min(U)),max(max(U))])
    drawnow
    pause(T)
end

end

