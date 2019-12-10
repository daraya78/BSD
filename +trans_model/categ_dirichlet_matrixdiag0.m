classdef categ_dirichlet_matrixdiag0 < handle
    properties
        prior = struct('categ_dirichlet',[]); %Sampling of Parameters
        posterior = struct('categ_dirichlet',[]); %Sampling of Parameters
        parsampl
        eelog    %Log Expectation
        expectation   %Expectation
        divkl    %Divergence
        ndim        
    end
    
    methods
        function self = categ_dirichlet_matrixdiag0(ndim)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            self.prior.categ_dirichlet.conc=[];
            self.posterior.categ_dirichlet.conc=[];
            self.eelog=[];
            self.divkl=[];
            self.expectation=[];
        end
        function self=parsamplfun(self,option)
            if option==1
                for j=1:self.ndim
                    aux=self.posterior.categ_dirichlet.conc(j,:);
                    aux(j,j)=0;
                    self.parsampl=self.posterior.categ_dirichlet.conc./repmat(sum(self.posterior.categ_dirichlet.conc,2),1,self.ndim);
                end
            elseif option==2
                for j=1:self.ndim
                    aux=self.prior.categ_dirichlet.conc(j,:);
                    aux(1,j)=0;
                    self.parsampl(j,:)=util.Dirichlet.dirichletRnd(aux);
                end
            elseif option==3
                for j=1:self.ndim
                    aux=self.posterior.categ_dirichlet.conc(j,:);
                    aux(j,j)=0;
                    self.parsampl=self.posterior.categ_dirichlet.conc./repmat(sum(self.posterior.categ_dirichlet.conc,2),1,self.ndim);
                end
            elseif option==4
                for j=1:self.ndim
                    aux=self.prior.categ_dirichlet.conc(j,:);
                    aux(j,j)=0;
                    self.parsampl=self.prior.categ_dirichlet.conc./repmat(sum(self.prior.categ_dirichlet.conc,2),1,self.ndim);
                end
            end
        end
        function re = sample(self,num)
            for k=1:num 
                for j=1:self.ndim
                    re(j,:,k)=util.Dirichlet.dirichletRnd(self.parsampl(j,:));
                end
            end
        end
        function re=prob(self,option,X)
        %No tiene sentido
       end
        function self=update(self,opttrain,X)
            if strcmp(opttrain,'EM')  %Used for E.M
                sxi=sum(X,3);
                self.parsampl=X./repmat(sum(X,2),[1 self.ndim]); 
            elseif strcmp(opttrain,'VB')   %Used for Variational Bayes
                self.posterior.categ_dirichlet.conc = sum(X,3) + self.prior.categ_dirichlet.conc;
                self.posterior.categ_dirichlet.conc = self.posterior.categ_dirichlet.conc.*(ones(self.ndim)-eye(self.ndim));
            end
        end
        function [p,plog]=expect(self,opttrain)
            if strcmp(opttrain,'EM')
                re=self.parsampl;
            elseif strcmp(opttrain,'VB')
                for j=1:self.ndim
                    self.eelog(j,:)=psi(self.posterior.categ_dirichlet.conc(j,:))-psi(sum(self.posterior.categ_dirichlet.conc(j,:))); % 10.66
                end
                plog=self.eelog;
                self.eelog=exp(self.eelog);
                p=self.eelog;
            end
        end
        function re=priorfull(self)
            re=1;
            if isempty(self.prior.categ_dirichlet.conc)
               re=0; 
            end
        end
        function re=priornoinf(self)
            self.prior.categ_dirichlet.conc=ones(self.ndim)-eye(self.ndim);
        end
        function expectfun(self)
            self.expectation=self.posterior.categ_dirichlet.conc./repmat(sum(self.posterior.categ_dirichlet.conc,2),1,self.ndim);
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
                categ_dirichlet.conck=[self.posterior.categ_dirichlet.conc(j,1:j-1) self.posterior.categ_dirichlet.conc(j,j+1:end)];
                categ_dirichlet.conc0=[self.prior.categ_dirichlet.conc(j,1:j-1) self.prior.categ_dirichlet.conc(j,j+1:end)];
                self.divkl=self.divkl+util.Dirichlet.KLDirich(categ_dirichlet.conck,categ_dirichlet.conc0);
            end
        end
    end
end
