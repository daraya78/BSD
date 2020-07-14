classdef mar3 < handle
   
    properties
        %prior = struct('coef_prec_gamma',[]);
        prior = struct('coef_normal',[],'prec_gamma',[],'coef_mask',[]);
        posterior = struct('coef_normal',[], 'prec_gamma',[]); 
        parsampl = struct('coef',[], 'prec',[]);
        expectation
        eelog       %Expectation
        divkl       %Divergence
        ndim        %Dimension
        nstates     %States Numbers
        order
    end
    
    methods
        function self=mar3(ndim,nstates,order)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            if exist('nstates','var'), self.nstates=nstates; else self.nstates=[]; end
            if exist('order','var'), self.order=order; else self.order=[]; end
            self.prior.coef_normal{1}.mean=[];
            self.prior.coef_normal{1}.prec_gamma.shape=[];
            self.prior.coef_normal{1}.prec_gamma.scale=[];
            self.prior.prec_gamma{1}.shape=[];
            self.prior.prec_gamma{1}.scale=[];
            self.posterior.coef_normal{1}.mean=[];
            self.posterior.coef_normal{1}.prec_gamma.shape=[];
            self.posterior.coef_normal{1}.prec_gamma.scale=[];
            self.posterior.coef_normal{1}.prec=[];
            self.posterior.prec_gamma{1}.shape=[];
            self.posterior.prec_gamma{1}.scale=[];
            self.eelog=[];
            self.divkl=[];
            self.expectation=[];
         end
        function self=parsamplfun(self,option)
            for j=1:self.nstates
                %PENDIENTE
                if option==1     % Posterior based
                %HOLD
                elseif option==2 % Prior based
                %HOLD
                elseif option==3 % Expectation Posterior based
                    %d=self.ndim;
                    %p=self.order;
                    self.parsampl.coef{j}=self.posterior.coef_normal{j}.mean;
                    shape=self.posterior.prec_gamma{j}.shape;
                    scale=self.posterior.prec_gamma{j}.scale;
                    self.parsampl.prec{j}=diag(shape.*scale);    
                end
            end
        end   
        function re = sample2(self,num,state,dataant)
            ord=self.order;
            if exist('dataant','var')
                x = dataant;
            else
                x =rand(ord,self.ndim);
            end
            w = self.parsampl.coef{state};
            noise=mvnrnd(zeros(1,self.ndim),inv(self.parsampl.prec{state}),num);
            A=util.spm.spm_unvec(w,zeros(self.ndim,self.ndim*ord));
            AT      = A';
            u      = [x; zeros(num,self.ndim)];
            if ord==1
                for k=2:num+1; 
                    x(1,:) = u(k-1,:)*AT;
                    u(k,:) = x + noise(k-1,:);
                end
            else
                for k=ord+1:num+ord; 
                    for j=1:ord;
                        x(j,:) = u(k-j,:)*AT((j-1)*self.ndim+1:j*self.ndim,:);
                    end
                    u(k,:) = sum(x)+noise(k-ord,:);
                    utemp(k,:)=sum(x);
                end
            end
            re=u(ord+1:end,:);
        end
        function [y noise] = sample(self,num,state,dataant)
            p=self.order;
            d=self.ndim;
            mask=self.prior.coef_mask;
            if exist('dataant','var')
                dataant=dataant(end-p+1:end,:);
            else
                dataant=rand(p,d);
            end
            if iscell(self.parsampl.coef{state})
                w0=util.convertermar(self.posterior.coef_normal{state}.mean,d,p,mask);
            else
               w0=self.parsampl.coef{state};
            end
            w = reshape(w0,d*p,d);
            noise=mvnrnd(zeros(1,self.ndim),inv(self.parsampl.prec{state}),num);
            y=dataant;
            for j=1:num
                dataaux=[y(end-p+1:end,:);ones(1,d)];
                x=util.regressor(dataaux,p);
                y=[y;x*w+noise(j,:)];
            end
            y(1:p,:)=[];
        end
        function [bnj logbnj]=prob(self,opttrain,y,state)
             if ~exist('state','var')
                k1=1;
                k2=self.nstates;
            else
                k1=state;
                k2=k1;
            end
            if strcmp(opttrain,'EM')       % Conventional (pdf) 
                bnj=mvnpdf(X,self.par{state}.mu,self.par{state}.Sigma);
                %Pendiente
            elseif strcmp(opttrain,'VB')  % Variational Bayes
                conta=1;
                [Xtot ytot]=util.regressor(y,self.order);
                N=size(ytot,1);
                d = self.ndim;
                p = self.order;
                mask=self.prior.coef_mask;
                maskm=reshape(mask,d*p,d)';
                for state=k1:k2
                    prob=-d/2*log(2*pi);
                    for j=1:d
                        wmean=self.posterior.coef_normal{state}.mean{j};
                        wprec=self.posterior.coef_normal{state}.prec{j};
                        yshape=self.posterior.prec_gamma{state}.shape(j);
                        yscale=self.posterior.prec_gamma{state}.scale(j);
                        maskw=maskm(j,:);
                        X=Xtot(:,find(maskw));
                        c2=1/2*(psi(yshape)+log(yscale));
                        aux1=-1/2*yshape*yscale;
                        aux2=(ytot(:,j)-X*wmean).^2+sum(X*inv(wprec).*X,2);
                        c3=aux1*aux2;
                        prob=prob+c2+c3;
                    end
                    logbnj(:,conta)=prob;
                    bnj(:,conta)=exp(prob);
                    conta=conta+1;
                end
            end
        end
        function update (self,opttrain,data,gammain,statereq)
            if ~exist('statereq','var')
                k1=1;
                k2=self.nstates;
            else
                k1=statereq;
                k2=k1;
            end
            if strcmp(opttrain,'EM') 

            elseif strcmp(opttrain,'VB') %Used for Variational Bayes
                for state=k1:k2
                    if ~exist('statereq','var')
                        gamma=gammain(:,state);
                    else
                        gamma=gammain;
                    end
                    if size(self.prior.coef_normal,2)>1
                        nprior=state;
                    else
                        nprior=1;
                    end
                    p=self.order;
                    d=self.ndim;
                    N=size(data,1);                % length of time series
                    Nef=sum(gamma);
                    [Xtot ytot]=util.regressor(data,p);
                    mask=self.prior.coef_mask;
                    maskm=reshape(mask,d*p,d)';
                    wshape0m=reshape(self.prior.coef_normal{nprior}.prec_gamma.shape,d*p,d)';
                    wscale0m=reshape(self.prior.coef_normal{nprior}.prec_gamma.scale,d*p,d)';
                    yshape0=self.prior.prec_gamma{nprior}.shape;
                    yscale0=self.prior.prec_gamma{nprior}.scale;
                    for j=1:d
                        y=ytot(:,j);
                        maskw=logical(maskm(j,:));
                        wshape0=wshape0m(j,maskw);
                        wscale0=wscale0m(j,maskw);
                        X=Xtot(:,maskw);
                        nr=sum(maskw);   
                        X=X.*repmat(sqrt(gamma),1,nr);
                        y=y.*sqrt(gamma);
                        if length(self.posterior.coef_normal)==self.nstates && (length(self.posterior.coef_normal{end}.mean)==d)
                            w_mean=self.posterior.coef_normal{state}.mean{j};
                            w_cov=inv(self.posterior.coef_normal{state}.prec{j});
                        else
                            % Initialition using ML (%Penny paper based)
                            w_ml{j}=pinv(X)*y;
                            %w_ml{j}=pinv(Xor).*repmat(gamma',nr,1)*y;
                            %y_pred=Xor*w_ml{j};
                            y_pred=X*w_ml{j};
                            e=y-y_pred;
                            %e=yor-y_pred;
                            %ycov_ml{j}=(e'*(e.*gamma))/Nef;  %Revisar si no es Nefectivo
                            
                            ycov_ml{j}=(e'*e)/Nef;
                            yprec_ml{j}=inv(ycov_ml{j});
                            w_covml{j}=kron(ycov_ml{j},inv(X'*X)); %Revisar
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            w_mean=w_ml{j};
                            w_cov=w_covml{j};
                        end
                        for i=1:1
                            % Update weight precisions
                            E_wj=0.5*w_mean'*w_mean;
                            wscalek=(wscale0.^-1+0.5*(w_mean'*w_mean+trace(w_cov))).^-1;  %Eqs 1
                            wshapek=wshape0+0.5;                                          %Eqs 1
                            mean_alpha=wshapek.*wscalek;
                            % Update noise precision posterior
                            ee=y-X*w_mean;
                            yshapek=yshape0(j)+Nef/2;                                       %Eqs 2
                            yscalek=(yscale0(j)^-1+0.5*(ee'*ee+sum(sum(X*w_cov.*X))))^-1; %Eqs 2
                            noise=yshapek*yscalek;
                            cov2(j)=noise;
                            % Update weight posterior
                            w_cov=inv(diag(mean_alpha)+noise*X'*X);                       %Eqs 3
                            w_prec=inv(w_cov);                                            %Eqs 3
                            w_mean=noise*w_cov*X'*y;
                        end
                        self.posterior.coef_normal{state}.mean{j}=w_mean;
                        self.posterior.coef_normal{state}.prec{j}=w_prec;
                        self.posterior.coef_normal{state}.prec_gamma.shape{j}=wshapek;
                        self.posterior.coef_normal{state}.prec_gamma.scale{j}=wscalek;
                        self.posterior.prec_gamma{state}.shape(j)=yshapek;
                        self.posterior.prec_gamma{state}.scale(j)=yscalek;
                    end
                end
            end
        end         
        function n=ndatafun(self,X)
            n=size(X,1)-self.order;
        end
        function priornoinf(self,type,X)
            if ~isempty(self.ndim)
                n=self.ndim;
                d=self.order;   
                if strcmp(type,'default')
                    self.prior.coef_normal{1}.prec_gamma.shape=0.001*ones(1,n*n*d);
                    self.prior.coef_normal{1}.prec_gamma.scale=1000*ones(1,n*n*d);
                    self.prior.coef_normal{1}.mean_normal.mean=zeros(1,n*n*d);
                    self.prior.coef_mask=ones(1,n*n*d);
                    self.prior.prec_gamma{1}.shape=0.001*ones(1,n*n*d);
                    self.prior.prec_gamma{1}.scale=1000*ones(1,n*n*d);
                elseif strcmp(type,'databased')
                    %HOLD
                end
            end
        end
        function cleanpos(self)
            self.posterior.coef_mask=[];
            self.posterior.prec_gamma=[];
            self.posterior.coef_normal=[];
            self.eelog=[];
        end
        function copy(self,o1)
            self.prior = o1.prior;
            self.posterior = o1.posterior;
            self.parsampl = o1.parsampl;
            self.eelog = o1.eelog;
            self.divkl = o1.divkl;
            self.ndim = o1.ndim;
            self.nstates = o1.nstates;
            self.order=o1.order;
        end
        function re=posteriorfull(self)
            re=0;
            if (length(self.posterior.coef_normal)==self.nstates) && (self.nstates~=0)
                re=1;
            end
        end
        function re=priorfull(self)
            re=1;
            coefsc=isempty(self.prior.coef_normal{1}.prec_gamma.scale);
            coefsh=isempty(self.prior.coef_normal{1}.prec_gamma.shape);    
            de=isempty(self.prior.prec_gamma{1}.shape);
            sc=isempty(self.prior.prec_gamma{1}.scale);
            if coefsc || coefsh || de || sc
               re=0; 
            end
        end
        function expectfun(self)
            for j=1:self.nstates
                self.expectation.prec{j}=self.posterior.prec_gamma{j}.scale.*self.posterior.prec_gamma{j}.shape;
                self.expectation.coef{j}=self.posterior.coef_normal{j};   
            end
        end
        function []=divklfun(self)
            D=0;
            d=self.ndim;
            p=self.order;
            mask=self.prior.coef_mask;
            maskm=reshape(mask,d*p,d)';
            Dyprec=0;
            Dwmean=0;
            Dwprec=0;
            for state=1:self.nstates
                Dgam=0;
                if size(self.prior.coef_normal,2)>1
                    nprior=state;
                else
                    nprior=1;
                end
                yshape0=self.prior.prec_gamma{nprior}.shape;
                yscale0=self.prior.prec_gamma{nprior}.scale;
                wmean02=reshape(self.prior.coef_normal{nprior}.mean_normal.mean,d*p,d);
                wshape02=reshape(self.prior.coef_normal{nprior}.prec_gamma.shape,d*p,d)';
                wscale02=reshape(self.prior.coef_normal{nprior}.prec_gamma.scale,d*p,d)';
                yshapek=self.posterior.prec_gamma{state}.shape;
                yscalek=self.posterior.prec_gamma{state}.scale;
                for j=1:d
                     maskw=logical(maskm(j,:));
                     nr=sum(maskw);
                     wshape0=wshape02(maskw);
                     wscale0=wscale02(maskw);
                     wmean0=wmean02(maskw);
                     wshapek=self.posterior.coef_normal{state}.prec_gamma.shape{j};
                     wscalek=self.posterior.coef_normal{state}.prec_gamma.scale{j};
                     wprec0=diag(wshapek.*wscalek);
                     wmeank=self.posterior.coef_normal{state}.mean{j};
                     wpreck=self.posterior.coef_normal{state}.prec{j};
                     Dwmean=Dwmean+util.Normal.spm_kl_normal(wmeank,inv(wpreck),wmean0,inv(wprec0));
                     Dyprec=Dyprec+util.Gamma.klgamma(yshapek(j),1/yscalek(j),yshape0(j),1/yscale0(j));
                     for k=1:nr
                        Dwprec=Dwprec+util.Gamma.klgamma(wshapek(k),1/wscalek(k),wshape0(k),1/wscale0(k));                         
                     end
                end
            end
            D=Dwmean+Dyprec+Dwprec;
            self.divkl=D;
        end
        function [flag]=trouble(self)
           flag=0;
           for state=1:self.nstates
                for j=1:size(self.posterior.coef_normal{1}.prec,2)
                    flag=flag+(rcond(self.posterior.coef_normal{state}.prec{j})<1e-10);
                end 
           end
           flag=flag>0;
        end
        function video=graf(self,opt,decodevar,video)
            if opt.graf==1
                fig=area(decodevar.gamma(1:2000,:));
                fig(1).Parent.Position=[0.13   0.117   0.775   0.81];
                drawnow
            end
            if opt.video==1
                   frame = getframe(gcf);
                   writeVideo(video,frame);
            end
        end
    end
end

