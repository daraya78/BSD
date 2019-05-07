function [decodevar] = init(self,X,opt)
    ndataefec=self.emis_model.ndatafun(X);
    if strcmp(opt.option,'kmeans')
        ncluster=self.nstates;
        stateseq=kmeans(X,ncluster,'replicates',opt.nrep,'MaxIter',opt.maxiter);
        %load stateseq
        [durationseq stateseqnorep]=util.durseq(stateseq');     
    else
        %%%%%%%%%%%%%RANDOM SAMPLE PARAMETER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        self.init(2);    
        %%%%%%%%%%GENERATE DATA AND DELETE PAR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [datatmp stateseq ndata stateseqnorep durationseq] = self.gen(ndataefec,'fixedndata','yes');
    end
    %%%%%%%%GENERATE GAMMA, XI AND ETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gamma=zeros(ndataefec,self.nstates);
    for j=1:self.nstates
        gamma(stateseq==j,j)=1;
    end
    if strcmp(class(self),'hsmm')
        nnorep=size(stateseqnorep,2);
        xi=zeros(self.nstates,self.nstates,nnorep);
        eta=zeros(opt.dmax,self.nstates,nnorep);
        for i=1:nnorep-1
            xi(stateseqnorep(i),stateseqnorep(i+1),i)=1;
            eta(durationseq(i),stateseqnorep(i),i)=1;
        end
        eta(durationseq(i+1),stateseqnorep(i+1),i+1)=1;
    else
       
        xi=zeros(self.nstates,self.nstates,ndataefec);
        for i=1:ndataefec-1
            xi(stateseq(i),stateseq(i+1),i)=1;
        end
        eta=[];
    end
    decodevar.gamma=gamma;
    decodevar.xi=xi;
    decodevar.eta=eta;
end