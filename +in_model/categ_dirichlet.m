classdef categ_dirichlet < handle
    properties
        prior = struct('categ_dirichlet',[]); %Sampling of Parameters
        posterior = struct('categ_dirichlet',[]); %Sampling of Parameters
        parsampl
        eelog    %Log Expectation
        expectation   %Expectation
        divkl       %Divergence
        ndim
        in_mar
    end
    
    methods
        function self = categ_dirichlet(ndim,in_mar)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            self.prior.categ_dirichlet.conc=[];
            self.posterior.categ_dirichlet.conc=[];
            if exist('in_mar','var'), self.in_mar=in_mar; else self.in_mar=[]; end
            self.eelog=[];
            self.divkl=[];
            self.expectation=[];
        end
        function self=parsamplfun(self,option)
           if option==1
                self.parsampl=util.Dirichlet.dirichletRnd(self.posterior.categ_dirichlet.conc);
           elseif option==2
                self.parsampl=util.Dirichlet.dirichletRnd(self.prior.categ_dirichlet.conc);
           elseif option==3
                self.parsampl=self.posterior.categ_dirichlet.conc/sum(self.posterior.categ_dirichlet.conc);
           elseif option==4
               self.parsampl=self.prior.categ_dirichlet.conc/sum(self.prior.categ_dirichlet.conc);
           end 
        end
        function re = sample(self,num)
            re=util.sample_discrete(self.parsampl,num);
        end 
        function [p,plog]=prob(self,opttrain,X)
            if strcmp(opttrain,'EM')      %Conventional
                %pi_0=self.pi_0; Pendiente
            elseif strcmp(opttrain,'VB')  %Variational Bayes
                plog=psi(self.posterior.categ_dirichlet.conc)-psi(sum(self.posterior.categ_dirichlet.conc,2));
                p=exp(p);
                %re=plog(X);
            end
        end
        function [p,plog]=expect(self,opttrain)
           if strcmp(opttrain,'EM')
               %p=self.parsampl; 
           elseif strcmp(opttrain,'VB')
               plog=psi(self.posterior.categ_dirichlet.conc)-psi(sum(self.posterior.categ_dirichlet.conc,2));
               self.eelog=exp(plog);
               p=self.eelog;
           end 
        end
        function self=update(self,opttrain,X)
           
            if strcmp(opttrain,'EM')
                %
            elseif strcmp(opttrain,'VB')
                %count=hist(X,1:self.ndim);
                %self.posterior.categ_dirichlet.conc = self.prior.categ_dirichlet.conc+count;      
                self.posterior.categ_dirichlet.conc = self.prior.categ_dirichlet.conc+X;      
                
            end
        end
        function re=priornoinf(self)
            self.prior.categ_dirichlet.conc=ones(1,self.ndim);
        end 
        function cleanpos(self)
            self.posterior.categ_dirichlet.conc=[];
            self.eelog=[];
        end
        function copy(self,i1)
            self.prior = i1.prior;
            self.posterior = i1.posterior;
            self.parsampl = i1.parsampl;
            self.eelog = i1.eelog;
            self.divkl = i1.divkl;
            self.ndim = i1.ndim;
        end  
        function re=priorfull(self)
            re=1;
            if isempty(self.posterior.categ_dirichlet.conc)
               re=0; 
            end
        end
        function expectfun(self)
            self.expectation=self.posterior.categ_dirichlet.conc/sum(self.posterior.categ_dirichlet.conc);
        end
        function re=divklfun(self)
            self.divkl=util.Dirichlet.KLDirich(self.posterior.categ_dirichlet.conc,self.prior.categ_dirichlet.conc);
        end
    end
end
