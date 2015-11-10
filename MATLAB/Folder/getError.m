function [ error ] = getError( U,u )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


k = size(U,2);
Err1 = zeros(1,k);
for i = 1:k
    Err1(i) = norm(u(:,i)-U(:,i) ,Inf)/max(abs(u(:,i)));
end

m = size(U,1);
Err2 = zeros(1,m);
for i = 1:m
    Err2(i) = norm(u(i,:)-U(i,:) ,Inf)/max(abs(u(i,:)));
end

error = min([Err1,Err2]); % rart?

end

