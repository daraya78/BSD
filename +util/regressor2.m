function [reg y]=regressor2(data,p)
N=size(data,1);
d=size(data,2);
reg=[];
for kp=p-1:-1:0
    for j=1:d
        aux=data(kp+1:N-(p-kp),j);
        reg=[reg aux];
    end
end
y=data([p+1:1:N],:);


