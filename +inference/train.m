function [out,sim_array, model, modelstruct modelstructarray] = train(self,data,opt_varargin)
    %Inputs
    % data:         Emission Data
    % opt_varargin: options change of default
    %Output
    % out:          Result Struct
    % ma:
    % sim_array:
    % model:
    n=size(opt_varargin,2);
    ent=0;
    for j=1:2:n
        if strcmp(opt_varargin{j},'opt')
            opt=opt_varargin{j+1};
            ent=1;
        end
    end
    if ent==0
    program_option; %Charge Default options
    end
    if opt.parallel
        numcores=feature('numcores');
        c=parcluster;
        c.NumWorkers=numcores*3;
    end
    if ~isempty(self.nstates) && self.nstates>0
        self.in_model.ndim=self.nstates;
        self.trans_model.ndim=self.nstates;
    end
    
    ma=[];
    for j=1:2:n
        opt=setfield(opt,opt_varargin{j},opt_varargin{j+1});
    end
    if ~self.emis_model.priorfull()
        if ~isstruct(data)
            self.emis_model.priornoinf(opt.prior,data);
        end
    end
    if strcmp(class(self),'hsmm')
        if ~self.dur_model.priorfull()
            self.dur_model.priornoinf(opt.prior,opt.dmax);
        end
    end
    if ~isempty(self.nstates) && (self.nstates~=0)
        if ~self.trans_model.priorfull()
            self.trans_model.priornoinf();
        end
        if ~self.in_model.priorfull()
            self.in_model.priornoinf();
        end
        if opt.parallel
            parfor j=1:opt.maxitersim
                model(j)=self.copytrain(0);
                [sim_array(j),modelstructarray(j)]=model(j).trainn(data,j,opt);
            end
        else
            for j=1:opt.maxitersim
                model(j)=self.copytrain(0);
                [sim_array(j),modelstructarray(j)]=model(j).trainn(data,j,opt);
            end
        end
        [a pos2]=max(cell2mat({sim_array.fe}));
        self.copy(model(pos2));
        out=sim_array(pos2);
        modelstruct=modelstructarray(pos2);
    else
        ma=-inf(1,opt.maxstates);
        for kstate=opt.minstates:opt.maxstates
            kstate
            self.nstates=kstate;
            if opt.parallel
                parfor j=1:opt.maxitersim
                    model(j,kstate)=self.copytrain(1,kstate);
                    if ~model(j,kstate).trans_model.priorfull()
                        model(j,kstate).trans_model.priornoinf();
                    end
                    if ~model(j,kstate).in_model.priorfull()
                        model(j,kstate).in_model.priornoinf();
                    end
                    [sim_array(j,kstate) modelstructarray(j,kstate)]=model(j,kstate).trainn(data,j,opt);
                end
            else
                for j=1:opt.maxitersim
                    model(j,kstate)=self.copytrain(1,kstate);
                    if ~model(j,kstate).trans_model.priorfull()
                        model(j,kstate).trans_model.priornoinf();
                    end
                    if ~model(j,kstate).in_model.priorfull()
                        model(j,kstate).in_model.priornoinf();
                    end 
                    [sim_array(j,kstate) modelstructarray(j,kstate)]=model(j,kstate).trainn(data,j,opt);
                end
            end
            array_ite(:,kstate)=cell2mat({sim_array(:,kstate).fe})';
            [ma(kstate) pos3(kstate)]=max(array_ite(:,kstate));
            if ma(kstate)<ma(kstate-1)
                break
            end
        end
        if ~(ma(kstate)<ma(kstate-1))
            self.copy(model(pos3(kstate),kstate));
            out=sim_array(pos3(kstate),kstate);
            modelstruct=modelstructarray(pos3(kstate),kstate);
        else
            self.copy(model(pos3(kstate-1),kstate-1));
            out=sim_array(pos3(kstate-1),kstate-1);
            modelstruct=modelstructarray(pos3(kstate-1),kstate-1);
        end
    end
    
    if ~isstruct(data)
        self.expectfun()
    else
        modelstruct.expectfun();
    end
    
    if strcmp(opt.seq,'maxgamma')
        if ~isstruct(data)
            [a pos2]=max(out.decodevar.gamma');
            out.stateseq=pos2;
            out.nstates=length(unique(pos2));
        else
            for kcond=1:size(out.decodevar.cond,2)
                for ksubj=1:size(out.decodevar.cond(kcond).subj,2)
                    for kblock=1:size(out.decodevar.cond(kcond).subj(ksubj).block,2)
                        gamma=out.decodevar.cond(kcond).subj(ksubj).block(kblock).decodevar.gamma;
                        [a pos2]=max(gamma');
                        out.stateseq.cond(kcond).subj(ksubj).block(kblock).stateseq=pos2;
                        out.nstates.cond(kcond).subj(ksubj).block(kblock).nstates=length(unique(pos2));
                    end
                end
            end
        end    
    elseif strcmp(opt.seq,'viterbi')
        if ~isstruct(data)
            out.stateseq=self.viterbi(data,opt);
            out.nstates=length(unique(out.stateseq));
        else
            out.stateseq=modelstruct.viterbi(data,opt);
            for kcond=1:size(out.decodevar.cond,2)
                for ksubj=1:size(out.decodevar.cond(kcond).subj,2)
                    for kblock=1:size(out.decodevar.cond(kcond).subj(ksubj).block,2)
                        nstates=length(unique(out.stateseq.cond(kcond).subj(ksubj).block(kblock).stateseq));
                        out.nstates.cond(kcond).subj(ksubj).block(kblock).nstates=nstates;
                    end
                end
            end
            
        end
    end
end      
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    