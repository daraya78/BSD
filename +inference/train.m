function [out,sim_array, model, modelstruct] = train(self,data,opt_varargin)
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
    ma=[];
    for j=1:2:n
        opt=setfield(opt,opt_varargin{j},opt_varargin{j+1});
    end
    if ~isempty(self.nstates) && (self.nstates~=0)
        for j=1:opt.maxitersim
            model(j)=self.copytrain(0);
            [sim_array(j),modelstruct]=model(j).trainn(data,j,opt);
        end
        [a pos2]=max(cell2mat({sim_array.fe}));
        self.copy(model(pos2));
        out=sim_array(pos2);
    else
        ma=-inf(1,opt.maxstates);
        for kstate=opt.minstates:opt.maxstates
            kstate
            for j=1:opt.maxitersim
                self.nstates=kstate;
                model(j,kstate)=self.copytrain(1,kstate);
                model(j,kstate).in_model.priornoinf();
                model(j,kstate).trans_model.priornoinf();
                [sim_array(j,kstate) modelstructma(j,kstate)]=model(j,kstate).trainn(data,j,opt);
                %if ~isempty(aux), auxmodelstructma(j,kstate)=aux;end
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
            modelstruct=modelstructma(pos3(kstate),kstate);
        else
            self.copy(model(pos3(kstate-1),kstate-1));
            out=sim_array(pos3(kstate-1),kstate-1);
            modelstruct=modelstructma(pos3(kstate-1),kstate-1);
        end
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
                        out.nstates.cond(kcond).subj(ksubj).block(kblock).stateseq=length(unique(pos2));
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
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    