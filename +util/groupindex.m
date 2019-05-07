function [re2 nre pos] = groupindex(cond,n,share)
    flag=zeros(n.ncond,max(n.nsubj),max(max(n.nblock)));
    for kcond=1:n.ncond
        for ksubj=1:n.nsubj(kcond)
            for kblock=1:n.nblock(kcond,ksubj)
                mcond=(1-share(1))*kcond+share(1);
                msubj=(1-share(2))*ksubj+share(2);
                mblock=(1-share(3))*kblock+share(3);
                if flag(mcond,msubj,mblock)==0
                    %re(mcond,msubj,mblock).data=[cond(kcond).subj(ksubj).block(kblock).data];
                    re(mcond,msubj,mblock,1,:)=[kcond ksubj kblock];
                    flag(kcond,ksubj,kblock)=1;
                else
                    m=size(re,4);
                    re(mcond,msubj,mblock,m+1,:)=[kcond ksubj kblock];
                end
                
            end
        end
    end
    conta=1;
    for kcond=1:size(re,1)
        for ksubj=1:size(re,2)
            for kblock=1:size(re,3)
                aux=squeeze(re(kcond,ksubj,kblock,:,:));
                if size(aux,2)==1,aux=aux';end
                if aux(1)~=0
                    re2{conta}=aux(aux(:,1)~=0,:);
                    conta=conta+1;
                end
            end
        end
    end
    nre=size(re2,2);
    
    
    for j=1:size(re2,2)
        ac=0;
        for k=1:size(re2{j},1)
            ac=ac+n.ndata(re2{j}(k,1),re2{j}(k,2),re2{j}(k,3));
            pos{j}(k)=ac;
        end
        
    end
    
    
end
   


