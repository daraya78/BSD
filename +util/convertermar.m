function [out1 out2]=convertermar(in,d,p,mask)
    if size(in,3)>1
        ma2=[];
        for j=1:p
            ma2=[ma2 in(:,:,j)];
            out2{j}=in(:,:,j);
        end
        out1=reshape(ma2',1,d*d*p); %Mascara
    elseif size(in,1)==1 & iscell(in)==0
        aux=reshape(in,d*p,d)';
        for j=1:p
            out1(:,:,j)=aux(:,(j-1)*d+1:j*d);
            out2{j}=out1(:,:,j);
        end
    elseif iscell(in)
        if nargin()==3 %Estructura de matlab
            ma2=[];
            for j=1:p
                out2(:,:,j)=in{j};
                ma2=[ma2 in{j}];
            end
            out1=reshape(ma2',1,d*d*p); %Mascara
        elseif nargin()==4  %estructura de hsmm
            out1=[];
            maskm=reshape(mask,d*p,d)';
            for j=1:d
                aux=zeros(1,d*p);
                aux(logical(maskm(j,:)))=reshape(in{j},1,sum(maskm(j,:)));
                out1=[out1 aux];
            end
            out2=[];
        end
            
    end
end     
        
            
            
        