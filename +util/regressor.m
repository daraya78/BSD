    function [reg y]=regressor(data,p)
                N=size(data,1);
                d=size(data,2);
                reg=[];
                for i=1:p,
                    tmpx=data(i:N-p+i,:);
                    reg=[reg,tmpx]; % Regresor Matrix
                end
                for i=1:p
                    start=(i-1)*d+1;
                    stop=start+d-1;
                    chunk(i).reg=reg(:,[start:1:stop]); %Se debe optimizar
                end
                reg=[];
                for i=p:-1:1,
                    reg=[reg chunk(i).reg];
                end
                % and remove last row
                Nrows=size(reg,1);
                reg=reg(1:Nrows-1,:);
                y=data([p+1:1:N],:);
        end
   