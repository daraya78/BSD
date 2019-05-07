function X=inmatr(Xma,X,ind,field,pos)
    
    if ~isobject(Xma)
        nind=size(ind,1);
        in=1;
        for j=1:nind
            pal=['X.cond(ind(' num2str(j) ',1)).subj(ind(' num2str(j) ',2)).block(ind(' num2str(j) ',3)).' field   '=Xma(in:pos(j),:);'];
            eval(pal);           
            in=pos(j)+1;
        end
    else
        nind=size(ind,1);
        in=1;
        for j=1:nind
            pal=['X(ind(' num2str(j) ',1),ind(' num2str(j) ',2),ind(' num2str(j) ',3)).' field '=Xma;'];
            eval(pal);
        end
    end
    
    
    
end