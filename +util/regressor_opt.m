    function [reg y]=regressor(data,p)
          
%     N=size(data,1);
%     d=size(data,2);
%     reg=[];
%     for i=1:p,
%         tmpx=data(i:N-p+i,:);
%         reg=[reg,tmpx]; % Regresor Matrix
%     end
% 
%     % the next two loops can be replaced using array operations
%     % using two lines of code (see bellow alternative code)
%     for i=1:p
%         start=(i-1)*d+1;
%         stop=start+d-1;
%         chunk(i).reg=reg(:,[start:1:stop]); %Se debe optimizar
%     end
%     reg=[];
%     for i=p:-1:1,
%         reg=[reg chunk(i).reg];
%     end
% 
% 
%     % Removing last row might not be necessary if the first
%     % loop is done as I suggested above. In the alternative
%     % code bellow I have taken that into account, but check
%     % this is true. If this is not true, then use the first
%     % loop as you have it now and include removing last row at
%     % the end of the alternative code bellow.
% 
%     % and remove last row
%     Nrows=size(reg,1);
%     reg=reg(1:Nrows-1,:);
%     y=data([p+1:1:N],:);


%% Code from Nelson reproducing previous code (same result, more efficient)
N=size(data,1);
d=size(data,2);
reg= zeros(N-p,d,p);
for i=0:1:p-1
    reg(:,:,i+1) = data(i+1:N-p+i,:); % Regresor Matrix
end

reg = reg(:,:,(p:-1:1));
reg = reshape(reg,N-p,d*p);
y=data([p+1:1:N],:);
                
        end
   