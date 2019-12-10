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
    for kcond=1:size(seq.nstates,1)
        for ksubj=1:size(seq.nstates.cond(kcond).subj,2)
            for kblock=1:size(seq.nstates.cond(kcond).subj(ksubj).block,2)
                nstates2=seq.nstates.cond(kcond).subj(ksubj).block(kblock).nstates;
                n=n+1;
                if nstates2>nstates
                    nstates=nstates2;
                end
            end
        end
    end
    
    f=figure
    conta=1;
    for kcond=1:size(seq.stateseq,1)
        for ksubj=1:size(seq.stateseq.cond(kcond).subj,2)
            for kblock=1:size(seq.stateseq.cond(kcond).subj(ksubj).block,2)
                subplot(2,ceil(n/2),conta);
                seq2=seq.stateseq.cond(kcond).subj(ksubj).block.stateseq;
                [a1 b11]=hist(seq2,1:nstates);
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


