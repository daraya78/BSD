classdef normal_normal_gamma < handle
   
    properties
        prior = struct('mean_normal',[], 'prec_gamma',[]); %Sampling of Parameters
        posterior = struct('mean_normal',[], 'prec_gamma',[]); %Sampling of Parameters
        parsampl = struct('mean',[], 'prec',[]); %Sampling of Parameters
        eelog       %log Expectation
        expectation      %Expectation
        divkl       %Divergence
        ndim        %Dimension
        nstates     %States Numbers
    end
    
    
    methods
        function self=normal_normal_gamma(ndim,nstates)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            if exist('nstates','var'), self.nstates=nstates; else self.nstates=[]; end
            self.prior.mean_normal{1}.mean=[];
            self.prior.mean_normal{1}.prec=[];
            self.prior.prec_gamma{1}.shape=[];
            self.prior.prec_gamma{1}.scale=[];
            self.posterior.mean_normal{1}.mean=[];
            self.posterior.mean_normal{1}.prec=[];
            self.posterior.prec_gamma{1}.shape=[];
            self.posterior.prec_gamma{1}.scale=[];
            self.eelog=[];
            self.expectation=[];
            self.divkl=[];
        end
        function self=parsamplfun(self,option)
            for j=1:self.nstates
                if option==1
                    self.parsampl.prec{j}=diag(gamrnd(self.posterior.prec_gamma{j}.shape,self.posterior.prec_gamma{j}.scale));
                    self.parsampl.mean{j}=mvnrnd(self.posterior.mean_normal{j}.mean,(self.posterior.mean_normal{j}.prec)^-1);
                elseif option==2
                    if size(self.prior.mean_normal,2)>1
                        self.parsampl.prec{j}=diag(gamrnd(self.prior.prec_gamma{j}.shape,self.prior.prec_gamma{j}.scale));
                        self.parsampl.mean{j}=mvnrnd(self.prior.mean_normal{j}.mean,(self.prior.mean_normal{j}.prec)^-1);
                    else
                        %flag=0;
                        %while (flag==0)
                            prec=diag(gamrnd(self.prior.prec_gamma{1}.shape,self.prior.prec_gamma{1}.scale));
                            flag=det(prec)>0.00001;
                        %end
                        self.parsampl.prec{j}=prec;
                        self.parsampl.mean{j}=mvnrnd(self.prior.mean_normal{1}.mean,(self.prior.mean_normal{1}.prec)^-1);
                    end
                elseif option==3
                    self.parsampl.prec{j}=diag(self.posterior.prec_gamma{j}.shape.*self.posterior.prec_gamma{j}.scale);
                    self.parsampl.mean{j}=self.posterior.mean_normal{j}.mean;
                elseif option==4
                    flag=size(self.prior.mean_normal,2)>1;
                    k=flag*j+1*(1-flag);
                    self.parsampl.prec{j}=diag(self.prior.prec_gamma{k}.shape.*self.prior.prec_gamma{k}.scale);
                    self.parsampl.mean{j}=self.prior.mean_normal{k}.mean;
                end
            end
        end
        function re = sample(self,num,state,dataant)
            re = mvnrnd(self.parsampl.mean{state},inv(self.parsampl.prec{state}),num);
        end
        function [p,plog]=prob(self,opttrain,X,state)
            if ~exist('state','var')
                k1=1;
                k2=self.nstates;
            else
                k1=state;
                k2=k1;
            end
            if strcmp(opttrain,'EM')       % Conventional (pdf)
                conta=1;
                for state=k1:k2
                    p(:,conta)=mvnpdf(X,self.parsampl.mean{state},inv(self.parsampl.prec{state}));
                    plog(:,conta)=log(p(:,conta));
                    conta=conta+1;
                end
            elseif strcmp(opttrain,'VB')  % Variational Bayes
                conta=1;
                for state=k1:k2
                    aux=sum(psi(self.posterior.prec_gamma{state}.shape)+log(self.posterior.prec_gamma{state}.scale));
                    XC = bsxfun(@minus, X, self.posterior.mean_normal{state}.mean);
                    E = trace( diag(self.posterior.prec_gamma{state}.shape.*self.posterior.prec_gamma{state}.scale)*(self.posterior.mean_normal{state}.prec)^-1)+...
                    sum((XC*diag(self.posterior.prec_gamma{state}.shape.*self.posterior.prec_gamma{state}.scale)).*XC,2);
                    plog(:,conta)=repmat(0.5*aux, size(X,1),1)- self.ndim*log(2*pi)- 0.5*E;
                    p(:,conta)=exp(plog(:,conta));
                    conta=conta+1;
                end
            end
        end
        function update (self,opttrain,X,gammain,statereq)
            if ~exist('statereq','var')
                k1=1;
                k2=self.nstates;
            else
                k1=statereq;
                k2=k1;
            end
            if strcmp(opttrain,'EM')    % Used for E.M
                for state=k1:k2
                     if ~exist('statereq','var')
                        gamma=gammain;
                    else
                        gamma=gammain(:,state);
                    end
                    scale=sum(gamma);
                    self.parsampl.mean{state}=sum(X.*repmat(gamma,[1 self.ndim]))/scale;
                    d=(X-ones(length(X),1)*self.parsampl.mean{state});
                    Sigma=(repmat(gamma,[1 self.ndim]).*d)'*d;  %Revisar
                    Sigma=Sigma/scale;
                    self.parsampl.prec{state}=inv(Sigma);
                end
            elseif strcmp(opttrain,'VB') %Used for Variational Bayes
                for state=k1:k2
                    if ~exist('statereq','var')
                        gamma=gammain(:,state);
                    else
                        gamma=gammain;
                    end
                    ndata=size(X,1);
                    Nk = sum(gamma,1);
                    Nk = Nk + 1e-10;
                    xbar = sum(bsxfun(@times, X, gamma)) / Nk;
                    if size(self.prior.mean_normal,2)>1
                        nprior=state;
                    else
                        nprior=1;
                    end
                    R0=self.prior.mean_normal{nprior}.prec;
                    m0=self.prior.mean_normal{nprior}.mean;
                    k0=self.prior.prec_gamma{nprior}.shape;
                    teta0=self.prior.prec_gamma{nprior}.scale;
                    
                    
                    if length(self.posterior.mean_normal)==self.nstates  
                        mk=self.posterior.mean_normal{state}.mean;
                        Rk=self.posterior.mean_normal{state}.prec;        
                    else
                        mk=sum(X.*repmat(gamma,1,self.ndim))./sum(gamma);
                        self.posterior.mean_normal{state}.mean=mk;
                        XC = bsxfun(@minus,X,xbar);
                        Sk = bsxfun(@times, XC, gamma)'*XC / Nk; %
                        Rk = Nk*inv(Sk);
                        self.posterior.mean_normal{state}.prec=Rk;
                    end
                    if Nk < 0.001 % extinguished
                        if size(self.prior.mean_normal,2)>1
                            self.posterior.mean_normal{state}.prec = self.prior.mean_normal{state}.prec; 
                            self.posterior.mean_normal{state}.mean = self.prior.mean_normal{state}.mean;
                            self.posterior.prec_gamma{state}.shape = self.prior.prec_gamma{state}.shape; 
                            self.posterior.prec_gamma{state}.scale = self.prior.prec_gamma{state}.scale; 
                        else
                            self.posterior.mean_normal{state}.prec = self.prior.mean_normal{1}.prec; 
                            self.posterior.mean_normal{state}.mean = self.prior.mean_normal{1}.mean;
                            self.posterior.prec_gamma{state}.shape = self.prior.prec_gamma{1}.shape; 
                            self.posterior.prec_gamma{state}.scale = self.prior.prec_gamma{1}.scale;
                        end
                    else
                        self.posterior.prec_gamma{state}.shape=Nk/2+k0;
                        aux=sum(repmat(gamma,1,self.ndim).*( (X-repmat(mk,ndata,1)).^2+repmat(1./diag(Rk)',ndata,1) ));
                        self.posterior.prec_gamma{state}.scale = 1./( 1/2*aux+1./teta0 );
                        kk=self.posterior.prec_gamma{state}.shape;
                        tetak=self.posterior.prec_gamma{state}.scale;
                        self.posterior.mean_normal{state}.prec=R0+Nk*diag(kk.*tetak);
                        self.posterior.mean_normal{state}.mean=(m0*R0+xbar*Nk*diag(kk.*tetak))*(self.posterior.mean_normal{state}.prec)^-1;
                    end               
                end
            end
        end 
        function n=ndatafun(self,X)
            n=size(X,1);
        end
        function priornoinf(self,type,X)
            n=self.ndim;
            if ~exist('type','var')
                type='default';
            end
            if strcmp(type,'default')
                self.prior.mean_normal{1}.mean=zeros(1,n);
                self.prior.mean_normal{1}.prec=diag(0.1*ones(1,n));
                self.prior.prec_gamma{1}.shape=0.001*ones(1,n);
                self.prior.prec_gamma{1}.scale=1000*ones(1,n);
                
            elseif strcmp(type,'databased')
                self.prior.mean_normal{1}.mean=mean(X);
                self.prior.mean_normal{1}.prec=diag(1./(((max(X)-min(X))/2).^2));
                self.prior.prec_gamma{1}.shape=0.001*ones(1,n);
                self.prior.prec_gamma{1}.scale=1000*ones(1,n); 
            end
        end
        function cleanpos(self)
            self.posterior.mean_normal=[];
            self.posterior.prec_gamma=[];
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
        end
        function re=posteriorfull(self)
            re=0;
            if (length(self.posterior.mean_normal)==self.nstates) && (self.nstates~=0)
                re=1;
            end
        end
        function re=priorfull(self)
            re=1;
            me=isempty(self.prior.mean_normal{1}.mean);
            pr=isempty(self.prior.mean_normal{1}.prec);
            sh=isempty(self.prior.prec_gamma{1}.shape);
            sc=isempty(self.prior.prec_gamma{1}.scale);
            if me || pr || sh || sc
               re=0; 
            end
        end
        function expectfun(self)
            for j=1:self.nstates
                self.expectation.prec{j}=diag(self.posterior.prec_gamma{j}.shape.*self.posterior.prec_gamma{j}.scale);
                self.expectation.mean{j}=self.posterior.mean_normal{j}.mean;
            end
        end
        function [flag]=trouble(self)
            flag=0;
        end
        function []=divklfun(self)
            D=0;
            for state=1:self.nstates
                 if size(self.prior.mean_normal,2)>1
                    R0=self.prior.mean_normal{state}.prec;
                    m0=self.prior.mean_normal{state}.mean;
                    k0=self.prior.prec_gamma{state}.shape;
                    teta0=self.prior.prec_gamma{state}.scale;
                else
                    R0=self.prior.mean_normal{1}.prec;
                    m0=self.prior.mean_normal{1}.mean;
                    k0=self.prior.prec_gamma{1}.shape;
                    teta0=self.prior.prec_gamma{1}.scale;
                end
                Rk=self.posterior.mean_normal{state}.prec;
                mk=self.posterior.mean_normal{state}.mean;
                kk=self.posterior.prec_gamma{state}.shape;
                tetak=self.posterior.prec_gamma{state}.scale; 
                for j=1:self.ndim         
                    D=D+util.Gamma.klgamma(kk(j),1/tetak(j),k0(j),1/teta0(j));   
                end                               
                D=D+util.Normal.spm_kl_normal(mk,inv(Rk),m0,inv(R0));  %Se debe revisar
            end
            self.divkl=D;
        end
    end
end

