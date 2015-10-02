function video(U,m,k,T)

for i = 1:k
    mesh(reshape(U(:,i),m,m),reshape(U(:,1),m,m)-max(max(U)))
    axis([0,m,0,m,min(min(U)),max(max(U))])
    %something colour
    %set(gcf,'Color',[1,0.4,0.6])
    drawnow
    pause(T)
end

end

