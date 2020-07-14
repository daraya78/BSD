function a1=ocupanciagraph(seq,tit)
    x0=100;
    y0=50;
    dx=200;
    dy=200;
if ~isstruct(seq)

    %f=figure;
    %f.Position=[x0 y0  x0+dx y0+dy];
    nstates=length(unique(seq));
    
    [a1 b11]=hist(seq,1:nstates);
    b=bar(b11,a1/length(seq)*100);
    x1=xlabel('BS','FontSize',10);
    y1=ylabel('BS Occupancy [%]','FontSize',10);
    y1.FontAngle='italic';
    x1.FontAngle='italic';
    grid
    title(tit)
    %axis([0 nstates+1 0 100]);

else
    
    nstates=0;
    n=0;
    uni=[];
    for kcond=1:size(seq.nstates,1)
        for ksubj=1:size(seq.nstates.cond(kcond).subj,2)
            for kblock=1:size(seq.nstates.cond(kcond).subj(ksubj).block,2)
                uni=[uni unique(seq.stateseq.cond(kcond).subj(ksubj).block(kblock).stateseq)];
            end
        end
    end
    n=max(unique(uni));
    f=figure
    conta=1;
    for kcond=1:size(seq.stateseq,1)
        for ksubj=1:size(seq.stateseq.cond(kcond).subj,2)
            for kblock=1:size(seq.stateseq.cond(kcond).subj(ksubj).block,2)
                subplot(2,ceil(5/2),conta);
                seq2=seq.stateseq.cond(kcond).subj(ksubj).block.stateseq;
                [a1 b11]=hist(seq2,1:n);
                b=bar(b11,a1/length(seq2)*100);
                bac(conta,:)=a1/length(seq2)*100;
                x1=xlabel('BS','FontSize',10);
                y1=ylabel('BS Occupancy [%]','FontSize',10);
                y1.FontAngle='italic';
                x1.FontAngle='italic';
                grid
                %title(tit)
                conta=conta+1;
                
            end
        end
    end
    f=figure;
    f.Position=[x0 y0  x0+dx y0+dy];
    b=bar(b11,sum(bac)/(conta-1));
    grid            
    x1=xlabel('BS','FontSize',10);
    y1=ylabel('BS Occupancy [%]','FontSize',10);
    y1.FontAngle='italic';
    x1.FontAngle='italic';


       
end

end


