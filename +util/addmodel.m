function [modelarray n iall] = addmodel(X,model)
    cond=X.cond;
    ncond=size(cond,2);
    n.ncond=(ncond);
    conta=1;
    for kcond=1:ncond
        nsubj=size(cond(kcond).subj,2);
        n.nsubj(kcond)=nsubj;
        for ksubj=1:nsubj
            nblock=size(cond(kcond).subj(ksubj).block,2);
            n.nblock(kcond,ksubj)=nblock;
            for kblock=1:nblock
                modelarray(kcond,ksubj,kblock)=model.copytrain(0);
                ndata=size(cond(kcond).subj(ksubj).block(kblock).data,1);
                n.ndata(kcond,ksubj,kblock)=ndata;
                iall(conta,:)=[kcond ksubj kblock];
                conta=conta+1;
            end
        end
    end
end

