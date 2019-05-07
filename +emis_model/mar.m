classdef mar < handle
   
    properties
        %prior = struct('coef_prec_gamma',[]);
        prior = struct('coef_prec_gamma',[],'prec_wishart',[]);
        posterior = struct('coef_mean_normal',[], 'prec_wishart',[]); 
        parsampl = struct('coef',[], 'prec',[]); 
        eelog       %Expectation
        divkl       %Divergence
        ndim        %Dimension
        nstates     %States Numbers
        order
    end
    
    methods
        
         function self=mar(ndim,nstates,order)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            if exist('nstates','var'), self.nstates=nstates; else self.nstates=[]; end
            if exist('order','var'), self.order=order; else self.order=[]; end
            self.priornoinf();
            self.posterior.coef_mean_normal{1}.mean=[];
            self.posterior.coef_mean_normal{1}.prec=[];
            self.posterior.coef_prec_gamma{1}.scale=[];
            self.posterior.coef_prec_gamma{1}.shape=[];
            self.posterior.prec_wishart{1}.degree=[];
            self.posterior.prec_wishart{1}.scale=[];
            self.eelog=[];
            self.divkl=[];
         end
        function self=parsamplfun(self,option)
            for j=1:self.nstates
                if option==1     % Posterior based

                elseif option==2 % Prior based
                    flag=size(self.prior.coef_prec_gamma,2)>1;
                    k=flag*j+1*(1-flag);
                    scale=self.prior.coef_prec_gamma{k}.scale;
                    shape=self.prior.coef_prec_gamma{k}.shape;
                    self.parsampl.alpha{j}=gamrnd(shape,scale);
                    self.parsampl.coef{j}=mvnrnd(zeros(1,self.ndim^2*self.order),inv(diag(self.parsampl.alpha{j})));
                    for i=1:10000
                        prec(:,:,i)=iwishrnd(eye(2),2);
                        lik(i)=abs(det(prec(:,:,i)))^-(self.ndim+1)/2;
                    end
                    pos=util.sample_discrete(lik,1);
                    self.parsampl.prec{j}=prec(:,:,pos);
                elseif option==3 % Expectation Posterior based
                    self.parsampl.coef{j}=self.posterior.coef_mean_normal{j}.mean;
                    degree=self.posterior.prec_wishart{j}.degree;
                    scale=self.posterior.prec_wishart{j}.scale;
                    self.parsampl.prec{j}=degree*scale; 
                end
            end
        end
            
        function re = sample(self,num,state,dataant)
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
                [x y]=util.regressor(y,self.order);
                for state=k1:k2
                    ndata=size(y,1);
                    d = self.ndim;
                    p = self.order;
                    mk = self.posterior.coef_mean_normal{state}.mean;
                    Rk = self.posterior.coef_mean_normal{state}.prec;
                    vk = self.posterior.prec_wishart{state}.degree;
                    Wk = self.posterior.prec_wishart{state}.scale;
                    alpha_scalek = self.posterior.coef_prec_gamma{state}.scale;
                    alpha_shapek = self.posterior.coef_prec_gamma{state}.shape;
                    loglambdaTilde = sum(psi(1/2*(vk + 1 -[1:d])))+d*log(2)+util.spm.spm_logdet(Wk);
                    w_ml=util.spm.spm_unvec(mk,zeros(d*p,d));
                    y_pred = x*w_ml;
                    e=y-y_pred;
                    invRk=inv(Rk);
                    invalpha=inv(diag(alpha_scalek.*alpha_shapek));
                    for j=1:ndata
                        ter1=trace(vk*Wk*e(j,:)'*e(j,:));
                        ter2=trace(vk*Wk*kron(eye(d),x(j,:))*invRk*kron(eye(d),x(j,:))');
                        logbnj(j,conta)=-d/2*log(2*pi)+0.5*loglambdaTilde-0.5*(ter1+ter2);
                    end
                    bnj(:,conta)=exp(logbnj(:,conta));
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
            
            if strcmp(opttrain,'EM') 

            elseif strcmp(opttrain,'VB') %Used for Variational Bayes
                for state=k1:k2
                    if ~exist('statereq','var')
                        gamma=gammain(:,state);
                    else
                        gamma=gammain;
                    end
                    if size(self.prior.coef_prec_gamma,2)>1
                        nprior=state;
                    else
                        nprior=1;
                    end
                    p=self.order;
                    d=self.ndim;
                    N=size(X,1);                % length of time series
                    [x y]=util.regressor(X,p);
                    x=x.*repmat(gamma,1,d*p);
                    y=y.*repmat(gamma,1,d);
                    Nef=sum(gamma);
                    k=p*d*d;
                    xp=pinv(-x);       %Pseudo inverse Moore-Penrose
                    inv_xtx=pinv(x'*x);
                    xtx=x'*x;
                    xty=x'*y;
                    % Get maximum likelihood solution
                    w_ml = xp*y;
                    % Swap signs to be consistent with paper (swap back at end !)
                    w_ml=-1*w_ml;
                    y_pred = x*w_ml;
                    e=y-y_pred;
                    noise_cov=(e'*e)/(Nef+p);
                    sigma_ml=kron(noise_cov,inv_xtx);
                    if length(self.posterior.coef_mean_normal)==self.nstates && ~isempty(self.posterior.coef_mean_normal{1}.mean)  
                        w_mean=reshape(self.posterior.coef_mean_normal{state}.mean,[p*d,d]);
                        w_cov=inv(self.posterior.coef_mean_normal{state}.prec);
                    else
                        w_mean=w_ml;
                        w_cov=sigma_ml;
                    end
                    max_iters=60;
                    w=zeros(p*d,d);
                    tol=0.0001;
                   %for it=1:max_iters,
                   it=1;
                        % Update weight precisions
                        mask=self.prior.coef_prec_mask{nprior}.mask;
                        groups=unique(mask);
                        shapekcom=zeros(1,k);
                        scalekcom=zeros(1,k);
                        for j=1:groups
                            maskgroup=(mask==j);
                            ind=find(maskgroup);
                            shape0=self.prior.coef_prec_gamma{nprior}.shape(ind(1));
                            scale0=self.prior.coef_prec_gamma{nprior}.scale(ind(1));
                            Ij=diag(maskgroup);
                            kj=sum(maskgroup);
                            E_wj=0.5*w_mean(:)'*Ij*w_mean(:); 
                            shapek=1 / (E_wj+0.5*trace(Ij*w_cov*Ij)+(1/shape0));
                            scalek=0.5*kj+scale0;
                            shapekcom(maskgroup)=shapek;
                            scalekcom(maskgroup)=scalek;
                        end
                        mean_alpha=shapekcom.*scalekcom;
                        % Update noise precision posterior
                        try
                        yy_pred=x*w_mean;
                        catch
                        keyboard
                        end
                        ee=y-yy_pred;
                        E_d_av=ee'*ee;
                        Omega = util.spm.spm_get_omega(p,d,w_cov,xtx);
                        E_d_av=E_d_av+Omega;
                        B=E_d_av+util.spm.spm_inv(self.prior.prec_wishart{nprior}.scale);
                        invB=util.spm.spm_inv(B);
                        a=Nef+p+self.prior.prec_wishart{nprior}.degree;
                        mean_lambda=a*invB;
                        % Update weight posterior
                        data_precision=kron(mean_lambda,xtx);
                        w_cov=util.spm.spm_inv(data_precision+diag(mean_alpha));
                        vec_w_mean=w_cov*data_precision*w_ml(:);
                        w_meanant=w_mean;
                        w_mean=reshape(vec_w_mean,p*d,d);
                        difw=sum(sum(abs(w_mean-w_meanant)));    
                        if (difw)< 10^-6
                            break
                        end    
                    %end
                    self.posterior.coef_mean_normal{state}.mean=util.spm.spm_vec(w_mean);
                    self.posterior.coef_mean_normal{state}.prec=inv(w_cov);
                    self.posterior.prec_wishart{state}.degree=a;
                    self.posterior.prec_wishart{state}.scale=invB;
                    self.posterior.coef_prec_gamma{state}.shape=shapekcom;
                    self.posterior.coef_prec_gamma{state}.scale=scalekcom;
                end
                for i=1:p,
                    start=(i-1)*d+1;
                    stop=(i-1)*d+1+(d-1);
                    % Transpose and swap signs for compatibility with spectral estimation function
                    mar.lag(i).a=w_mean(start:stop,:)';
                end
            end
            
        end
        
        function n=ndatafun(self,X)
            n=size(X,1)-self.order;
        end
        
        function priornoinf(self,X,varargin)
            
            if ~isempty(self.ndim)
                opt.prior='default';   
                for j=1:2:(nargin-2)
                 opt=setfield(opt,varargin{j},varargin{j+1});
                end
                n=self.ndim;
                d=self.order;
                if strcmp(opt.prior,'default')
                    self.prior.coef_prec_gamma{1}.shape=1000*ones(1,n*n*d);
                    self.prior.coef_prec_gamma{1}.scale=0.001*ones(1,n*n*d);
                    self.prior.coef_prec_mask{1}.mask=ones(1,n*n*d);
                    %self.prior.prec_wishart{1}.degree=n;
                    %self.prior.prec_wishart{1}.scale=0.001*eye(n);
                    self.prior.prec_wishart{1}.degree=self.ndim;
                    self.prior.prec_wishart{1}.scale=eye(n);

                
                elseif strcmp(opt.prior,'databased')
                    %HOLD
                end
            end
        end
        
        
        function cleanpos(self)
            self.posterior.coef_mean_normal=[];
            self.posterior.prec_wishart=[];
            self.posterior.coef_prec_gamma=[];
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
            if (length(self.posterior.coef_mean_normal)==self.nstates) && (self.nstates~=0)
                re=1;
            end
        end
        
        function []=divklfun(self)
            D=0;
            
            for state=1:self.nstates
                Dgam=0;
                if size(self.prior.coef_prec_gamma,2)>1
                    coefshape0=self.prior.coef_prec_gamma{state}.shape;
                    coefscale0=self.prior.coef_prec_gamma{state}.scale;
                    degree0=self.prior.prec_wishart{state}.degree;
                    scale0=self.prior.prec_wishart{state}.scale;
                else
                    coefshape0=self.prior.coef_prec_gamma{1}.shape;
                    coefscale0=self.prior.coef_prec_gamma{1}.scale;
                    degree0=self.prior.prec_wishart{1}.degree;
                    scale0=self.prior.prec_wishart{1}.scale;
                end
                coefprec0=diag(coefshape0.*coefscale0);
                coefmeank=self.posterior.coef_mean_normal{state}.mean;
                coefpreck=self.posterior.coef_mean_normal{state}.prec;
                coefshapek=self.posterior.coef_prec_gamma{state}.shape;
                coefscalek=self.posterior.coef_prec_gamma{state}.scale;
                degreek=self.posterior.prec_wishart{state}.degree;
                scalek=self.posterior.prec_wishart{state}.scale;
                
                Dnor=util.spm.spm_kl_eig_normal(coefmeank,inv(coefpreck),inv(coefprec0));
                Dwis=util.Wishart.spm_kl_wishart(degreek,inv(scalek),degree0,inv(scale0));
                for j=1:self.order*self.order*self.ndim
                    Dgam=Dgam+util.Gamma.klgamma(coefscalek(j),1/coefshapek(j),coefscale0(j),1/coefshape0(j));
                end
                D=D+Dnor+Dgam+Dwis;
            end
            self.divkl=D;
        end
        
        function graf(self,optiongraf,data,elib,cycle,stateseq,alphak)
            if optiongraf==1
                for j=1:self.nstates
                    Prec(:,:,j)=diag(self.hpark{j}.teta.*self.hpark{j}.k);    
                    %Prec(:,:,j)=self.obs_model.hpark{j}.v*self.obs_model.hpark{j}.W;
                    m(j,:)=self.hpark{j}.m;
                    R(:,:,j)=self.hpark{j}.R;     
                end
                %   util.plotFn(self.data,self.trans_model.alphak,m,Prec,nor,cycle,self.stateseq);
                util.plotFn(data(:,100:101),alphak,m(:,100:101),Prec(100:101,100:101,:),elib,cycle,stateseq);
                %util.plotFn(self.data(:,1:2),self.trans_model.alphak,m(:,1:2),R(1:2,1:2,:),elibre(cycle).value,cycle,self.stateseq);  
            end
        end
        
        
    end
end

