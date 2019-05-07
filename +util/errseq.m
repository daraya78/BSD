function [err1,err2,dif]=errseq(data,seqorig,m1,m2)
    dec1=m1.decode(data);
    [a pos1]=max(dec1.gamma');
    %pos1=m1.viterbi(data);
    seq1=util.renameseq(seqorig,pos1);
    if strcmp(class(m1.emis_model),'emis_model.mar')==1
        n=m1.emis_model.order+1;
    else
        n=1;
    end
    err1=(sum(seq1~=seqorig(n:end))/length(seqorig(n:end)))*100;
    if exist('m2','var')
        dec2=m2.decode(data);
        [a pos2]=max(dec2.gamma');
        seq2=util.renameseq(seqorig,pos2);
        err2=(sum(seq2~=seqorig(n:end))/length(seqorig(n:end)))*100;
        dif=err2-err1;
    end
    
end


    
