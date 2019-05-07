function [seq compara ] = renameseq(orig,estim)
orig_un=unique(orig);
estim_un=unique(estim);
n1=length(orig_un);
n2=length(estim_un);
seq=zeros(1,n1);
final=10;
compara=[];
while ~isempty(estim_un)
    ranking=[];
    for j=1:size(estim_un,2)
        maskj=(estim==estim_un(j));
        [counts val]=hist(orig(maskj),0:1:10);
        counts(1)=[];val(1)=[];
        [ma pos]=max(counts);
        ranking(j,1)=val(pos);
        ranking(j,2)=ma;
    end
    [ma pos]=max(ranking(:,2));
    if ranking(pos,1)==0
        seq(estim==estim_un(pos))=final;
        compara=[compara;[final estim_un(pos)]];
        estim_un(pos)=[];
        final=final-1;
        
        disp('Revisar')
    else
        seq(estim==estim_un(pos))=ranking(pos,1);
        compara=[compara;[ranking(pos,1) estim_un(pos)]];
        estim_un(pos)=[];
        orig(orig==ranking(pos,1))=0;
        
    end    
    

end
    
[aux1 pos3]=sort(compara(:,1));
compara=compara(pos3,:);   
    
    
   
end



