function [stateseq] = initvb(self,X,opt,repi)
    ndataefec=self.emis_model.ndatafun(X);
    if strcmp(opt.initoption,'kmeans')
        ncluster=self.nstates;
        stateseq=util.kmeans2015(X,ncluster,'replicates',opt.nrep,'MaxIter',opt.maxiter);
        stateseq=stateseq';
        %load stateseq
        [durationseq stateseqnorep]=util.durseq(stateseq);     
    elseif strcmp(opt.initoption,'random')
        %%%%%%%%%%%%%RANDOM SAMPLE PARAMETER%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %self.init(2);    
        %%%%%%%%%%GENERATE DATA AND DELETE PAR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %self.trans_model.parsampl=[0.8 0.1 0.1;0.1 0.8 0.1;0.1 0.1 0.8]; %Provisorio
        %for j=1:100
        %    [datatmp stateseq ndata stateseqnorep durationseq] = self.gen(ndataefec,'fixedndata','yes');
        %    if length(unique(stateseq))==self.nstates,break;end
        %end
        stateseq=[];
        while(length(unique(stateseq))~=self.nstates)
            stateseq=util.unirndseq(ndataefec,self.nstates,opt.dmax);
        end
        %load sqaux;
        %stateseq=sqaux;
    elseif strcmp(opt.initoption,'random2')
        stateseq=[];
        while(length(unique(stateseq))~=self.nstates)
            stateseq=util.unirndseq(ndataefec,self.nstates,opt.dmax);
        end

        
    elseif strcmp(opt.initoption,'seqprior')
        seqs=opt.seqs;
        nseqs=size(seqs,1);
        if nseqs==1
            if size(size(seqs),2)==3
                stateseq(:,:)=squeeze(seqs(1,:,self.nstates))';
            else
                stateseq=seqs';
            end
        else
            for j=1:nseqs
                opt2=opt;
                model(j)=self.copytrain(0);
                opt2.maxcyc=3;
                opt2.initoption='seqprior';
                opt2.seqs=seqs(j,:,self.nstates);
                %opt2.verbose='no';
                outdata(j)=model(j).trainn(X,opt2);
            end
            [a b]=max([outdata.fe]);
            stateseq=seqs(b,:,self.nstates)';
        end
    elseif strcmp(opt.initoption,'seqprior_mar')
        %if size(opt.seqs,3)==1
        %    stateseq=squeeze(opt.seqs(:,self.nstates))';
        %else
        %    stateseq=squeeze(opt.seqs(repi,:,self.nstates));
        %end
        stateseq=opt.seqs{self.nstates}(repi,:);
    end
     
    %%%%%%%%GENERATE GAMMA, XI AND ETA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %decodevar=util.decodevarinit(stateseq',class(self),opt.dmax,opt);
end