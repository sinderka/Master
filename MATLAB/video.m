function video(U,m,k,T)

for i = 3:k
    mesh(reshape(U(:,i),m,m))
    axis([0,m,0,m,-1,1])
    %something colour
    drawnow
    pause(T)
end

end

