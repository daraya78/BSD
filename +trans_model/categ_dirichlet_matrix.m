classdef categ_dirichlet_matrix < handle
    properties
        prior = struct('categ_dirichlet',[]); %Sampling of Parameters
        posterior = struct('categ_dirichlet',[]); %Sampling of Parameters
        parsampl
        eelog    %Expectation
        divkl       %Divergence
        ndim        
    end
    
    methods
        function self = categ_dirichlet_matrix(ndim)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            self.eelog=[];
            self.divkl=[];
            if ~isempty(ndim),self.priornoinf();end
        end
        function self=parsamplfun(self,option)
            if option==1
                for j=1:self.ndim
                    self.parsampl(j,:)=util.Dirichlet.dirichletRnd(self.posterior.categ_dirichlet.conc(j,:));
                end
            elseif option==2
                for j=1:self.ndim
                    self.parsampl(j,:)=util.Dirichlet.dirichletRnd(self.prior.categ_dirichlet.conc(j,:));
                end
            elseif option==3
                    self.parsampl=self.posterior.categ_dirichlet.conc./repmat(sum(self.posterior.categ_dirichlet.conc,2),1,self.ndim);
            elseif option==4
                    self.parsampl=self.prior.categ_dirichlet.conc./repmat(sum(self.prior.categ_dirichlet.conc,2),1,self.ndim);     
            end
        end
        function re = sample(self,num)
            for k=1:num 
                for j=1:self.ndim
                    re(j,:,k)=util.Dirichlet.dirichletRnd(self.parsampl(j,:));
                end
            end
        end
        function re=prob(self,opttrain,X)
        %No tiene sentido
       end
        function self=update(self,opttrain,X)
            if strcmp(opttrain,'EM')  %Used for E.M
                sxi=sum(X,3);
                self.parsampl=X./repmat(sum(X,2),[1 self.ndim]); 
            elseif strcmp(opttrain,'VB')   %Used for Variational Bayes
                self.posterior.categ_dirichlet.conc = sum(X,3) + self.prior.categ_dirichlet.conc;
            end
        end
        function [p,plog]=expect(self,opttrain)
            if strcmp(opttrain,'EM')
                p=self.parsampl;
            elseif strcmp(opttrain,'VB')
                for j=1:self.ndim
                    self.eelog(j,:)=psi(self.posterior.categ_dirichlet.conc(j,:))-psi(sum(self.posterior.categ_dirichlet.conc(j,:))); % 10.66
                end
                plog=self.eelog;
                self.eelog=exp(self.eelog);
            end
            p=self.eelog;
        end
        function re=priornoinf(self)
            self.prior.categ_dirichlet.conc=ones(self.ndim);
        end
        function cleanpos(self)
            self.posterior.categ_dirichlet.conc=[];
            self.eelog=[];
        end
        function copy(self,t1)
            self.prior = t1.prior;
            self.posterior = t1.posterior;
            self.parsampl = t1.parsampl;
            self.eelog = t1.eelog;
            self.divkl = t1.divkl;
            self.ndim = t1.ndim;
        end 
        function divklfun(self)
            self.divkl=0;
            for j=1:self.ndim
                categ_dirichlet.conck=self.posterior.categ_dirichlet.conc(j,:);
                categ_dirichlet.conc0=self.prior.categ_dirichlet.conc(j,:);
                self.divkl=self.divkl+util.Dirichlet.KLDirich(categ_dirichlet.conck,categ_dirichlet.conc0);
            end
        end
    end
end
