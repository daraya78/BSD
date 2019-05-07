function [redur reseq]=calculadur(seq)
    n=size(seq,2);
    redur=[];
    reseq=[];
    
    
    dur=1;
    for j=2:n
        if seq(j)==seq(j-1)
            dur=dur+1;
            if j==n
                redur=[redur dur];
                reseq=[reseq seq(j)];
            end
        else
            redur=[redur dur];
            dur=1;
            reseq=[reseq seq(j-1)];
            if j==n
                redur=[redur 1];
                reseq=[reseq seq(j)];
            end
        end
    end
end

            