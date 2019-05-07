function re=outmatr(Xma,X,ind,field,ncat)
    re=[];
    nind=size(ind,1);
    for j=1:nind
        pal=['X.cond(ind(' num2str(j) ',1)).subj(ind(' num2str(j) ',2)).block(ind(' num2str(j) ',3)).' field];
        if ncat~=4
            re=cat(ncat,re,eval(pal));          
        else
            aux=eval(pal);
            re=cat(1,re,aux(1,:));
        end
    end
end