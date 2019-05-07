function [decodevar] = decodevarinit(stateseq,type,opt)
        dmax=opt.dmax;
        ndataefec=size(stateseq,1);
        nstates=size(unique(stateseq),1);
        gamma=zeros(ndataefec,nstates);
        eta=[];
        durcount=zeros(nstates,dmax);
        for j=1:nstates
            gamma(stateseq==j,j)=1;
        end
        if strcmp(type,'hsmm')
            %[durationseq,stateseqnorep]=util.calculadur(stateseq);
            %nnorep=size(stateseqnorep,2);
            %xi=zeros(nstates,nstates,nnorep);
            %eta=zeros(dmax,nstates,nnorep);
            
            xi=zeros(nstates,nstates,ndataefec);
            if opt.deceta==1
                eta=zeros(dmax,nstates,ndataefec);
            end
            for i=1:ndataefec-1
                if (stateseq(i)~=stateseq(i+1))
                    xi(stateseq(i),stateseq(i+1),i)=1;
                    d=1;
                    if i==ndataefec-1
                        if opt.deceta==1
                            eta(1,stateseq(i+1),i)=1;
                        end
                        durcount(stateseq(i+1),1)=durcount(stateseq(i+1),1)+1;
                    else
                        while(stateseq(i+1)==stateseq(i+1+d)&& (i+1+d)<(ndataefec))
                            d=d+1;
                        end
                        d=min(d,dmax);
                        if (i+1+d)==ndataefec
                            if opt.deceta==1
                                eta(d+1,stateseq(i+1),i)=1;
                            end
                            durcount(stateseq(i+1),d+1)=durcount(stateseq(i+1),d+1)+1;
                        else
                            if opt.deceta==1
                                eta(d,stateseq(i+1),i)=1;
                            end
                            durcount(stateseq(i+1),d)=durcount(stateseq(i+1),d)+1;
                        end
                    end
                end
            end
            if opt.deceta==1
                eta=eta(1:dmax,:,:);
            end
        elseif strcmp(type,'hmm')   
            xi=zeros(nstates,nstates,ndataefec);
            for i=1:ndataefec-1
                xi(stateseq(i),stateseq(i+1),i)=1;
            end
        end
        decodevar.gamma=gamma;
        decodevar.xi=xi;
        decodevar.eta=eta;
        decodevar.durcount=durcount';
end

