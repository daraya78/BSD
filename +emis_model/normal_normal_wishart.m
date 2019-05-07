classdef normal_normal_wishart < handle
   
    properties
        prior = struct('mean_normal',[], 'prec_wishart',[]); %Sampling of Parameters
        posterior = struct('mean_normal',[], 'prec_wishart',[]); %Sampling of Parameters
        parsampl = struct('mean',[], 'prec',[]); %Sampling of Parameters
        eelog       %Expectation
        divkl       %Divergence
        ndim        %Dimension
        nstates     %States Numbers
    end
    
    methods
        
         function self=normal_normal_wishart(ndim,nstates)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            if exist('nstates','var'), self.nstates=nstates; else self.nstates=[]; end
            self.priornoinf();
            self.posterior.mean_normal{1}.mean=[];
            self.posterior.mean_normal{1}.prec=[];
            self.posterior.prec_wishart{1}.degree=[];
            self.posterior.prec_wishart{1}.scale=[];
            self.eelog=[];
            self.divkl=[];
         end

        function self=parsamplfun(self,option)
            for j=1:self.nstates
                if option==1
                    self.parsampl.prec{j}=iwishrnd(self.posterior.prec_wishart{j}.scale,self.posterior.prec_wishart{j}.degree);
                    self.parsampl.mean{j}=mvnrnd(self.posterior.mean_normal{j}.mean,(self.posterior.mean_normal{j}.prec)^-1);
                elseif option==2
                    if size(self.prior.mean_normal,2)>1
                        self.parsampl.prec{j}=iwishrnd(self.prior.mean_normal{j}.scale,self.prior.mean_normal{j}.degree);
                        self.parsampl.mean{j}=mvnrnd(self.prior.mean_normal{j}.mean,(self.prior.mean_normal{j}.prec)^-1);
                    else
                        flag=0;
                        while (flag==0)
                            prec=iwishrnd(self.prior.prec_wishart{1}.scale,self.prior.prec_wishart{1}.degree);
                            flag=det(prec)>0.00001;
                        end
                        self.parsampl.prec{j}=prec;
                        self.parsampl.mean{j}=mvnrnd(self.prior.mean_normal{1}.mean,(self.prior.mean_normal{1}.prec)^-1);
                    end
                elseif option==3
                    self.parsampl.prec{j}=self.posterior.prec_wishart{j}.scale*self.posterior.prec_wishart{j}.degree;
                    self.parsampl.mean{j}=self.posterior.mean_normal{j}.mean;   
                elseif option==4
                    flag=size(self.prior.mean_normal,2)>1;
                    k=flag*j+1*(1-flag);
                    self.parsampl.prec{j}=self.prior.prec_wishart{k}.scale*self.prior.prec_wishart{k}.degree;
                    self.parsampl.mean{j}=self.prior.mean_normal{k}.mean;  
                end
            end
        end
         
            
        function re = sample(self,num,state)
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
                bnj=mvnpdf(X,self.par{state}.mu,self.par{state}.Sigma);
                %Pendiente
            elseif strcmp(opttrain,'VB')  % Variational Bayes 
                conta=1;
                n=self.ndim;
                for state=k1:k2
                    degree=self.posterior.prec_wishart{state}.degree;
                    scale=self.posterior.prec_wishart{state}.scale;
                    m=self.posterior.mean_normal{state}.mean;
                    prec=self.posterior.mean_normal{state}.prec;
                    aux=sum(psi(1/2*(degree + 1 - [1:n]))) + ...
                    self.ndim*log(2) + util.spm.spm_logdet(scale);
                    XC = bsxfun(@minus, X, m);
                    E = trace( (degree*scale)*(prec)^-1)+degree*sum((XC*scale).*XC,2);
                    plog(:,conta)=repmat(0.5*aux, size(X,1),1)- n*log(2*pi)- 0.5*E;
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
            if strcmp(opttrain,'EM') % Used for E.M
                %scale=sum(gamma);
                %self.par{state}.mu=sum(X.*repmat(gamma,[1 self.ndim]))/scale;
                %d=(X-ones(length(X),1)*self.par{state}.mu);
                %self.par{state}.Sigma=(repmat(gamma,[1 self.ndim]).*d)'*d;
                %self.par{state}.Sigma=self.par{state}.Sigma/scale;    %Hold div n-1
            elseif strcmp(opttrain,'VB') %Used for Variational Bayes   
                for state=k1:k2
                    if ~exist('statereq','var')
                        gamma=gammain(:,state);
                    else
                        gamma=gammain;
                    end
                    ndata=size(X,1);
                    Nk = sum(gamma,1); % 10.51
                    Nk = Nk + 1e-10;
                    xbar=sum(bsxfun(@times, X, gamma)) / Nk;
                    if size(self.prior.mean_normal,2)>1
                        nprior=state;
                    else
                        nprior=1;
                    end
                
                    R0=self.prior.mean_normal{nprior}.prec;
                    m0=self.prior.mean_normal{nprior}.mean;
                    v0=self.prior.prec_wishart{nprior}.degree;
                    W0=self.prior.prec_wishart{nprior}.scale;
               
                    if length(self.posterior.mean_normal)==self.nstates  
                        mk=self.posterior.mean_normal{state}.mean;
                        Rk=self.posterior.mean_normal{state}.mean;
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
                            self.posterior.prec_wishart{state}.degree = self.prior.prec_wishart{state}.degree; 
                            self.posterior.prec_wishart{state}.scale = self.prior.prec_wishart{state}.scale; 
                        else
                            self.posterior.mean_normal{state}.prec = self.prior.mean_normal{1}.prec; 
                            self.posterior.mean_normal{state}.mean = self.prior.mean_normal{1}.mean;
                            self.posterior.prec_wishart{state}.degree = self.prior.prec_wishart{1}.degree; 
                            self.posterior.prec_wishart{state}.scale = self.prior.prec_wishart{1}.scale;                          
                        end
                    else               
                        XC = bsxfun(@minus,X,xbar);
                        Sk = bsxfun(@times, XC, gamma)'*XC / Nk;
                        self.posterior.prec_wishart{state}.degree=Nk+v0;
                        %Rk=self.posterior.mean_normal{1}.prec;
                        %self.posterior.prec_wishart{state}.scale=inv((W0^-1)+Nk*Sk+Nk*trace(inv(Rk)));
                        self.posterior.prec_wishart{state}.scale=inv((W0^-1)+Nk*Sk);                        
                        vk=self.posterior.prec_wishart{state}.degree;
                        Wk=self.posterior.prec_wishart{state}.scale;
                        self.posterior.mean_normal{state}.prec=R0+Nk*(vk*Wk);
                        self.posterior.mean_normal{state}.mean=(m0*R0+xbar*Nk*(vk*Wk))*(self.posterior.mean_normal{state}.prec)^-1;
                    end
                end
                    
            end
        end
        function n=ndatafun(self,X)
            n=size(X,1);
        end
        function priornoinf(self,X,varargin)
            opt.prior='default';   
            for j=1:2:(nargin-2)
                opt=setfield(opt,varargin{j},varargin{j+1});
            end
            n=self.ndim;
            if strcmp(opt.prior,'default')
                self.prior.mean_normal{1}.mean=zeros(1,n);
                %%%Provisorio
                if isempty(n)
                    n=0;
                end
                self.prior.mean_normal{1}.prec=0.01*eye(n);
                self.prior.prec_wishart{1}.degree=n;
                self.prior.prec_wishart{1}.scale=0.001*eye(n);
            elseif strcmp(opt.prior,'databased')
                self.prior.mean_normal{1}.mean=mean(X);
                self.prior.mean_normal{1}.prec=inv(cov(X)/self.ndim)
                self.prior.prec_wishart{1}.scale=inv(cov(X)/self.ndim);
                self.prior.prec_wishart{1}.degree=self.ndim;   
            end
        end
        function cleanpos(self)
            self.posterior.mean_normal=[];
            self.posterior.prec_wishart=[];
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
        function []=divklfun(self)
            D=0;
            for state=1:self.nstates
                if size(self.prior.mean_normal,2)>1    
                    R0=self.prior.mean_normal{state}.prec;
                    m0=self.prior.mean_normal{state}.mean;
                    v0=self.prior.prec_wishart{state}.degree;
                    W0=self.prior.prec_wishart{state}.scale;
                else
                    R0=self.prior.mean_normal{1}.prec;
                    m0=self.prior.mean_normal{1}.mean;
                    v0=self.prior.prec_wishart{1}.degree;
                    W0=self.prior.prec_wishart{1}.scale;
                end
                
                Rk=self.posterior.mean_normal{state}.prec;
                mk=self.posterior.mean_normal{state}.mean;
                vk=self.posterior.prec_wishart{state}.degree;
                Wk=self.posterior.prec_wishart{state}.scale;
                D=D+util.Normal.spm_kl_normal(mk,inv(Rk),m0,inv(R0))+util.Wishart.spm_kl_wishart(vk,inv(Wk),v0,inv(W0));
            end
            self.divkl=D;
        end
    end
end

