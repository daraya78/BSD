function topographstruct(outhsmm,data3,chanlocs)
    
    
    nstates=0;
    n=0;
    for kcond=1:size(outhsmm.nstates,1)
        for ksubj=1:size(outhsmm.nstates.cond(kcond).subj,2)
            for kblock=1:size(outhsmm.nstates.cond(kcond).subj(ksubj).block,2)
                nstates2=outhsmm.nstates.cond(kcond).subj(ksubj).block(kblock).nstates;
                n=n+1;
                if nstates2>nstates
                    nstates=nstates2;
                end
            end
        end
    end
    
    seq=[];
    for kcond=1:size(outhsmm.stateseq,1)
        for ksubj=1:size(outhsmm.stateseq.cond(kcond).subj,2)
            for kblock=1:size(outhsmm.stateseq.cond(kcond).subj(ksubj).block,2)
                seq=[seq outhsmm.stateseq.cond(kcond).subj(ksubj).block(kblock).stateseq];
            end
        end
    end
    
    
    
%Y =zscore(data3.X);
Y =data3.X;

%Y=data3.X-mean(data3.X);
X=zeros(length(seq),nstates);
for k=1:nstates
    X(seq==k,k) = 1;
end
%X = zscore(X);

[B, W] = util.spm.spm_robust_glm(Y, X, 1,[]);

f=figure;
f.Position=[27         344        1245         232];
ran=[min(min(B)) max(max(B))];
ran=[0.8 2.65]
%ran=[0.8 0.5];

for j=1:nstates
    subplot(2,nstates/2,j)
    %ran=[min(B(j,:)) max(B(j,:))];
    topoplot(B(j,:), chanlocs, 'electrodes','off','maplimits',ran);
    %topoplot(B(j,:), chanlocs, 'electrodes','off');
    
    colorbar()
end
ran
end