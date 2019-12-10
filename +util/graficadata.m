function graficadata(dataseg)
nstates=0;
n=0;
ma1=0;
ma2=0;
mi1=10000;
mi2=10000;
nsubj=size(dataseg.cond(1).subj,2);
for ksubj=1:nsubj
    kk=dataseg.cond(1).subj(ksubj).block(1).data;
    if max(kk(:,1))>ma1
        ma1=max(kk(:,1));
    end
    if min(kk(:,1))<mi1
        mi1=min(kk(:,1));
    end
    if max(kk(:,2))>ma2
        ma2=max(kk(:,2));
    end
    if min(kk(:,2))<mi2
        mi2=min(kk(:,2));
    end
end

figure
for j=1:5
        subplot(2,ceil(nsubj/2),j)
        kk=dataseg.cond(1).subj(j).block(1).data;
        plot(kk(:,1),kk(:,2),'x')
        hold on
        axis([mi1 ma1 mi2 ma2])
        %axis([-5 10 -6 6])
    end
end
