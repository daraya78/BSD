classdef categ_dirichlet < handle
   
    properties
        prior = struct('categ_dirichlet',[]);
        posterior = struct('categ_dirichlet',[]);
        parsampl = struct('conc',[]);
        eelog       %Expectation
        divkl       %Divergence
        ndim        %Dimension
        nstates     %States Numbers
    end
    
    methods
        function self=categ_dirichlet(ndim,nstates)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            if exist('nstates','var'), self.nstates=nstates; else self.nstates=[]; end
            self.priornoinf();
            self.posterior.categ_dirichlet{1}.conc=[];
            self.divkl=[];
        end
        function self=parsamplfun(self,option)
          for j=1:self.nstates
                if option==1
                    self.parsampl.conc{j}=util.Dirichlet.dirichletRnd(self.posterior.categ_dirichlet{j}.conc);
                elseif option==2
                    flag=size(self.prior.categ_dirichlet,2)>1;
                    k=flag*j+1*(1-flag);  
                    self.parsampl.conc{j}=util.Dirichlet.dirichletRnd(self.prior.categ_dirichlet{k}.conc);
                elseif option==3
                    self.parsampl.conc{j}=self.posterior.categ_dirichlet{j}.conc/sum(self.posterior.categ_dirichlet{j}.conc);
                elseif option==4
                    flag=size(self.prior.mean_normal,2)>1;
                    k=flag*j+1*(1-flag);                    
                    self.parsampl.conc{j}=self.prior.categ_dirichlet{j}.conc/sum(self.prior.categ_dirichlet{j}.conc);
                end 
          end
        end
        function re = sample(self,num,state)
            re=util.sample_discrete(self.parsampl.conc{state},num);
        end
        function [p2,plog2]=prob(self,opttrain,X,state)
            if ~exist('state','var')
                k1=1;
                k2=self.nstates;
            else
                k1=state;
                k2=k1;
            end
            if strcmp(opttrain,'EM')       % Conventional (pdf)

            elseif strcmp(opttrain,'VB')  % Variational Bayes
                conta=1;
                for state=k1:k2
                    plog=psi(self.posterior.categ_dirichlet{state}.conc)-psi(sum(self.posterior.categ_dirichlet{state}.conc,2));
                    p=exp(plog);
                    p2(:,conta)=p;
                    plog2(:,conta)=plog;
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
               
            elseif strcmp(opttrain,'VB') %Used for Variational Bayes
                for state=k1:k2
                    if ~exist('statereq','var')
                        gamma=gammain(:,state);
                    else
                        gamma=gammain;
                    end
                    if size(self.prior.categ_dirichlet,2)>1
                        nprior=state;
                    else
                        nprior=1;
                    end
                    self.posterior.categ_dirichlet{state}.conc = self.prior.categ_dirichlet{nprior}.conc+gamma';
                end
            end
        end
        function priornoinf(self,X,varargin)
            opt.prior='default';   
            for j=1:2:(nargin-2)
                opt=setfield(opt,varargin{j},varargin{j+1});
            end
            n=self.ndim;
            if strcmp(opt.prior,'default')
                self.prior.categ_dirichlet{1}.conc=ones(1,200) ;
            elseif strcmp(opt.prior,'databased')
            end
        end
        function cleanpos(self)
            self.posterior.categ_dirichlet=[];
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
        function []=divklfun(self)
            D=0;
            for state=1:self.nstates
                 if size(self.prior.categ_dirichlet,2)>1
                    conc0=self.prior.categ_dirichlet{state}.conc;
                else
                    conc0=self.prior.categ_dirichlet{1}.conc;
                end
                conck=self.posterior.categ_dirichlet{1}.conc;
                D=D+util.Dirichlet.KLDirich(conck,conc0);
            end
            self.divkl=D;
        end
    end
end

