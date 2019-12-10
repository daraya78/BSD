classdef multimodel < handle 
    properties
        matrixmodel
        n
        iall
        iemis
        itrans
        iin
        idur
    end
    methods
        function self = multimodel(matrixmodel,n,iall,iemis,itrans,iin,idur)
            if exist('matrixmodel','var'), self.matrixmodel=matrixmodel; else self.matrixmodel=[]; end
            if exist('n','var'), self.n=n; else self.n=[]; end
            if exist('iall','var'), self.iall=iall; else self.iall=[]; end
            if exist('iemis','var'), self.iemis=iemis; else self.iemis=[]; end
            if exist('itrans','var'), self.itrans=itrans; else self.itrans=[]; end
            if exist('iin','var'), self.iin=iin; else self.iin=[]; end
            if exist('idur','var'), self.idur=idur; else self.idur=[]; end
        end
       
        function [decstruct] = decode(self,X,opt)  
             for j=1:size(self.iall,1)
                model=self.matrixmodel(self.iall(j,1),self.iall(j,2),self.iall(j,3));
                Xaux=X.cond(self.iall(j,1)).subj(self.iall(j,2)).block(self.iall(j,3)).data;
                %auxdec=inference.hsmmresidual_decodelog(self,X,varargin)
                if ~exist('opt','var')
                    auxdec=model.decode(Xaux);
                else
                    auxdec=model.decode(Xaux,opt);
                end
                decstruct.cond(self.iall(j,1)).subj(self.iall(j,2)).block(self.iall(j,3)).decodevar=auxdec;
             end
        end
        function [decstruct] = viterbi(self,X,varargin)
             for j=1:size(self.iall,1)
                model=self.matrixmodel(self.iall(j,1),self.iall(j,2),self.iall(j,3));
                Xaux=X.cond(self.iall(j,1)).subj(self.iall(j,2)).block(self.iall(j,3)).data;
                [sq delta] = model.viterbi(Xaux);
                decstruct.cond(self.iall(j,1)).subj(self.iall(j,2)).block(self.iall(j,3)).stateseq=sq;
             end
        end        
        function init(self,option)
             for j=1:size(self.iall,1)
                model=self.matrixmodel(self.iall(j,1),self.iall(j,2),self.iall(j,3));
                model.init(option);
             end
        end   
        function expectfun(self)
            for j=1:size(self.iall,1)
                model=self.matrixmodel(self.iall(j,1),self.iall(j,2),self.iall(j,3));
                model.expectfun();
            end
        end
        function [out, sim_array, hsmm2, modelstruct]=train(self,data,varargin)
            [out, sim_array, hsmm2, modelstruct]=inference.train(self,data,varargin);
        end
        function [out modelstruct]=trainn(self,data,varargin)
            [out modelstruct]=inference.trainn(self,data,varargin);
        end
        %function cleanpos(self)
        %end
        
        %function copy(self,hsmm1)
        %end
        
        %function hsmm2=copytrain(self,optcleanpos,nstates)
        %end

        function output=elibrefun(self,decstruct)
            divklemis=0;
            for j=1:size(self.iemis,2)
                kc=self.iemis{j}(1,1);
                ks=self.iemis{j}(1,2);
                kb=self.iemis{j}(1,3);
                self.matrixmodel(kc,ks,kb).emis_model.divklfun();
                divklemis=divklemis+self.matrixmodel(kc,ks,kb).emis_model.divkl;
            end
            divkltrans=0;
            for j=1:size(self.itrans,2)
                kc=self.itrans{j}(1,1);
                ks=self.itrans{j}(1,2);
                kb=self.itrans{j}(1,3);
                self.matrixmodel(kc,ks,kb).trans_model.divklfun();
                divkltrans=divkltrans+self.matrixmodel(kc,ks,kb).trans_model.divkl;
            end
            divkldur=0;
            if strcmp(class(self.matrixmodel),'hsmm')
                for j=1:size(self.idur,2)
                    kc=self.idur{j}(1,1);
                    ks=self.idur{j}(1,2);
                    kb=self.idur{j}(1,3);
                    self.matrixmodel(kc,ks,kb).dur_model.divklfun();
                    divkldur=divkldur+self.matrixmodel(kc,ks,kb).dur_model.divkl;
                end
            end
            divklin=0;
            for j=1:size(self.iin,2)
                kc=self.iin{j}(1,1);
                ks=self.iin{j}(1,2);
                kb=self.iin{j}(1,3);
                self.matrixmodel(kc,ks,kb).in_model.divklfun();
                divklin=divklin+self.matrixmodel(kc,ks,kb).in_model.divkl;
            end
            loglik=0;
            for j=1:size(self.iall,1)
                kc=self.iall(j,1);
                ks=self.iall(j,2);
                kb=self.iall(j,3);
                loglik=loglik+decstruct.cond(kc).subj(ks).block(kb).decodevar.loglik;
            end
            output.divklemis=divklemis;
            output.divkltrans=divkltrans;
            output.divkldur=divkldur;
            output.divklin=divklin;
            output.loglik=loglik;
            output.totaldivkl=-divkltrans-divklemis-divkldur-divklin;
            output.fe=loglik+output.totaldivkl;
        end 
        
        
        
        
        
        
    end
end

