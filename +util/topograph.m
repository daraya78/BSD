function topograph(outhsmm,data,chanlocs)
    if ~exist('data','var')
        ruta1='P:\BrainDynamicsLab\Software\HSMM\Sim_paper\Realdata_13Nov\';
        file11='s7b1 clean broadband ENV.mat';
        file22='s7b1 clean broadband.mat';
        load([ruta1 file11]);
        load([ruta1 file22]);
        data=ADATA;
        chanlocs=EEG.chanlocs;
    end
    data5=data';
    seq=outhsmm.stateseq;
    nstates=length(unique(seq));
    upoints=1:10000;
    data6=data5(upoints,:);
    Y = zscore(data6(upoints,:));
    seqest = seq(upoints);
    Xmat= zeros(length(seqest),nstates);
    for j=1:nstates
        Xmat(seqest==j,j) = 1;
    end
        [ beta1, sigma, E, V ] = mvregress(Xmat,Y,'algorithm','cwls');

    f=figure;
    if nstates<5
        ny=1;nx=nstates;
    elseif nstates<9
        ny=2;nx=ceil(nstates/2);
    elseif nstates<13
        ny=3;nx=ceil(nstates/3);
    elseif nstates<21
        ny=4;nx=ceil(nstates/4);
    end  
    
    for j=1:nstates
        subplot(ny,nx,j)
        topoplot(beta1(j,:), chanlocs, 'electrodes','off');
    end

end