function [out modelarray] = trainn(self,X,repi,opt_varargin)
    if ~exist('repi','var'),repi=[];end
    n=0;
    modelarray=nan;
    if isempty(opt_varargin)
        program_option;
    elseif isstruct(opt_varargin{1})
        opt=opt_varargin{1};
       
    else
        aux_func.defaultopt;
        n=size(opt_varargin,2);
    end
    for j=1:2:n
        opt=setfield(opt,opt_varargin{j},opt_varargin{j+1});
    end
    
    if strcmp(opt.prepro,'meanvarnor')
        ndata=size(X,1);
        X=(X-repmat(mean(X),ndata,1))./repmat(sqrt(var(X)),ndata,1); 
    end
    
    if ~isstruct(X)
        %INITIAL CONDITION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%INIT E-STEP RANDOM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if self.emis_model.posteriorfull();
            decodevar=self.decode(X,opt);
        else 
            seqaux = util.initvb(self,X,opt,repi);
            %load stateseq
            %seqaux=stateseq;
            decodevar=util.decodevarinit(seqaux',class(self),opt);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dif=inf;
        cycle=0;
        
        %%%provisorio%%%%%%%%%%%%%%%
        %load decodevarbueno
        %decodevar=decodevar2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %save decodevarbueno decodevar
        if opt.video==1
                video = VideoWriter([opt.file '.avi']);
                video.FrameRate=2;
                open(video);
        end
        
        
        while dif>opt.tol && cycle<opt.maxcyc
            cycle=cycle+1;     
            %save anterior
            %area(decodevar.gamma)
            %drawnow
            %%%%%%%%%%%%%%%%%%%%%%M-STEP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Upgrade of Observation Model
            self.emis_model.update(opt.train,X,decodevar.gamma);
            if self.emis_model.trouble()
                freearray(cycle).fe=-inf;
                break;
            end
            if strcmp(class(self),'hsmm')
                %durcount=sum(decodevar.eta,3);
                %self.dur_model.update(2,(1:opt.dmax)',durcount);
                self.dur_model.update(opt.train,(1:opt.dmax)',decodevar.durcount);
            end
            self.trans_model.update(opt.train,decodevar.xi); %Upgrade of Transition Model 
            self.in_model.update(opt.train,decodevar.gamma(1,:)); %Upgrade of initial Model
            %%%%%%%%%%%%%%%%%%%%%%%E-STEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            decodevar=self.decode(X,opt);
            
            %save(['decodevar' num2str(cycle) '.mat'],'decodevar')
            %load(['decodevar' num2str(cycle) '.mat'])
            
            %subplot(2,1,1)
            %a1=area(decodevar.gamma(1:800,:));
            %title(['Iteración: ' num2str(cycle)])
            %a1(1).FaceColor='r';a1(2).FaceColor='g';a1(3).FaceColor='b';
            %subplot(2,1,2)
            %a2=area(exp(decodevar.Emi(:,1:800))'./sum(exp(decodevar.Emi(:,1:800)'),2));
            %a2(1).FaceColor='r';a2(2).FaceColor='g';a2(3).FaceColor='b';
            %drawnow
            
            
            %%%%%%%%%%%%%%%%%%FREE ENERGY ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%
            freearray(cycle)=self.elibrefun(decodevar);
            %%%%%%%%%%%%%%%%%%VERBOSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(opt.verbose,'yes')
                if cycle==1
                    fprintf('\nRepetition:%2i - Iteration Start \n',repi);
                end
                fprintf('It:%3i \tFree Energy:%3.5f \n',cycle,freearray(cycle).fe);
            end
            %%%%%%%%%%%%%%%%%%%FREE ENERGY CONDITION CHECK %%%%%%%%%%%%%%%%%%%%
            if cycle>1
                dif=(freearray(cycle).fe-freearray(cycle-1).fe);
                if dif< 0 %-1e-5
                    %disp('ERROR')ce
                    %save posterior
                    %dif
                    dif=inf;
                    %keyboard
                end                
            end
            if opt.video==1
                video=self.emis_model.graf(opt,decodevar,video);
            else
                if opt.graf>0
                    self.emis_model.graf(opt,decodevar);
                end
            end
        end
        

        out.fe=freearray(end).fe;
        out.loglik=freearray(end).loglik;
        out.totaldivkl=freearray(end).totaldivkl;
        out.fedetail=freearray;
        decodevar.xi=[]; %heavy size
        out.decodevar=decodevar;
        out.niter=cycle;
        if opt.video==1
            close(video)
        end
        
        
    else
        %%%%%%%%%%%%%%%%%%%INIT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        if ~self.emis_model.priorfull()
            [modelarr n iall] = util.addmodel(X,self);
            X2=util.outmatr([],X,iall,'data',1);
            self.emis_model.priornoinf(opt.prior,X2);
        end
        [modelarr n iall] = util.addmodel(X,self);
        [iemis nemis poe]=util.groupindex(X,n,opt.shareemis);
        [itrans ntrans pot]=util.groupindex(X,n,opt.sharetrans);
        [iin nin poin]=util.groupindex(X,n,opt.sharein);
        [idur ndur podur]=util.groupindex(X,n,opt.sharedur);
        modelarray=multimodel(modelarr,n,iall,iemis,itrans,iin,idur);
        seqstruct=[];
        for j=1:nemis
            Xemi{j}=util.outmatr([],X,iemis{j},'data',1);
            model=modelarray.matrixmodel(iemis{j}(1,1),iemis{j}(1,2),iemis{j}(1,3));
            flag=1;
            while(flag>0)
                seqemi = util.initvb(model,Xemi{j},opt);
                ini=1;
                flag=0;
                for k=1:length(poe{j})
                   flag=flag+(length(unique(seqemi(ini:poe{j}(k))))~=self.nstates); 
                   ini=poe{j}(k);
                end
            end
                seqstruct=util.inmatr(seqemi',seqstruct,iemis{j},'seq',poe{j});
        end
        for j=1:size(iall,1)
            seqaux=seqstruct.cond(iall(j,1)).subj(iall(j,2)).block(iall(j,3)).seq;
            auxdec=util.decodevarinit(seqaux,class(self),opt);
            decstruct.cond(iall(j,1)).subj(iall(j,2)).block(iall(j,3)).decodevar=auxdec; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dif=inf;
        cycle=0;
        while dif>opt.tol && cycle<opt.maxcyc
            cycle=cycle+1;
            for j=1:nemis
                model=modelarray.matrixmodel(iemis{j}(1,1),iemis{j}(1,2),iemis{j}(1,3));
                gammatmp=util.outmatr([],decstruct,iemis{j},'decodevar.gamma',1);
                model.emis_model.update(opt.train,Xemi{j},gammatmp);
                modelarray.matrixmodel=util.inmatr(model.emis_model,modelarray.matrixmodel,iemis{j},'emis_model'); 
            end

            for j=1:ntrans
                model=modelarray.matrixmodel(itrans{j}(1,1),itrans{j}(1,2),itrans{j}(1,3));
                xitmp=util.outmatr([],decstruct,itrans{j},'decodevar.xi',3);
                model.trans_model.update(opt.train,xitmp);
                modelarray.matrixmodel=util.inmatr(model.trans_model,modelarray.matrixmodel,itrans{j},'trans_model');
            end
            
            for j=1:nin
                model=modelarray.matrixmodel(iin{j}(1,1),iin{j}(1,2),iin{j}(1,3));
                intmp=util.outmatr([],decstruct,iin{j},'decodevar.gamma',4);
                intmp=sum(intmp,1);
                model.in_model.update(opt.train,intmp);
                modelarray.matrixmodel=util.inmatr(model.in_model,modelarray.matrixmodel,iin{j},'in_model');
            end
                 
            if strcmp(class(self),'hsmm')
                for j=1:ndur
                    model=modelarray.matrixmodel(idur{j}(1,1),idur{j}(1,2),idur{j}(1,3));
                    durcounttmp=util.outmatr([],decstruct,idur{j},'decodevar.durcount',3);
                    durcounttmp=sum(durcounttmp,3);
                    model.dur_model.update(opt.train,(1:opt.dmax)',durcounttmp);
                    modelarray.matrixmodel=util.inmatr(model.dur_model,modelarray.matrixmodel,idur{j},'dur_model');
                end
            end
            %%%%STEP E %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            decstruct=modelarray.decode(X,opt);         
            %%%%%%%%%%%%%%%%%%FREE ENERGY ESTIMATION %%%%%%%%%%%%%%%%%%%%%%
            totaldivkl=0;
            totalfe=0;
            loglik=0;
            freearray(cycle)=modelarray.elibrefun(decstruct);

            %%%%%%%%%%%%%%%%%%VERBOSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if strcmp(opt.verbose,'yes')
                if cycle==1
                    fprintf('\nIteration Start \n');
                end
                fprintf('It:%3i \tFree Energy:%3.5f \n',cycle,freearray(cycle).fe);
            end
            %%%%%%%%%%%%%%%%%%%FREE ENERGY CONDITION CHECK %%%%%%%%%%%%%%%%%%%%
            if cycle>1
                dif=(freearray(cycle).fe-freearray(cycle-1).fe);
                if dif< 0 %-1e-5
                    dif=inf;
                end
            end
        end
        out.fe=freearray(end).fe;
        out.loglik=freearray(end).loglik;
        out.totaldivkl=freearray(end).totaldivkl;
        out.fedetail=freearray;
        out.decodevar=decstruct;
        out.niter=cycle;
        
    end
    
    
    
    
    
    
    
    
    
    
    
    