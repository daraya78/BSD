classdef hmm < handle
    
    properties
        ndim
        nstates
        in_model
        trans_model
        emis_model
    end

    methods
        function self = hmm(ndim,nstates,emis_model2,in_model2,trans_model2)
            if exist('ndim','var'), self.ndim=ndim; else self.ndim=[]; end
            if exist('nstates','var'), self.nstates=nstates; else self.nstates=[]; end
            if exist('emis_model2','var')
                self.emis_model=emis_model2;
            else
                e=emis_model.mar(self.ndim,self.nstates,2);
                %e=emis_model.normal_normal_gamma(self.ndim,self.nstates);
                self.emis_model=e;
                self.emis_model.priornoinf();
            end
            if exist('in_model2','var')
                self.in_model=in_model2;
            else
                i=in_model.categ_dirichlet(self.nstates);
                self.in_model=i;
                if ~isempty(self.nstates),self.in_model.priornoinf();end
            end
            if exist('trans_model2','var')
                self.trans_model=trans_model2;
            else
                t=trans_model.categ_dirichlet_matrix(self.nstates);
                self.trans_model=t;
                if ~isempty(self.nstates),self.trans_model.priornoinf();end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%% HMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [decodevar] = decode(self,X,varargin)
            [decodevar] = inference.hmm_decodelog(self,X,varargin);
        end
        function [currentState, logP] = viterbi(self,data,varargin)
            [currentState, logP] = inference.hmmviterbi(self,data,varargin);
        end
        function [data,seq,ndata,reseq,redur] = gen(self,ndata,varargin)   
            opt=[];
            for j=1:2:(nargin-2)
                opt=setfield(opt,varargin{j},varargin{j+1});
            end
            nextstate_distr = self.in_model.parsampl;
            data = [];
            for j=1:ndata
                seq(j) = util.sample_discrete(nextstate_distr);
                data=[data; self.emis_model.sample(1,seq(j))];
                nextstate_distr = self.trans_model.parsampl(seq(j),:)';
            end
            ndata=size(data,1);
            [redur,reseq]=util.durseq(seq);
        end
        function init(self,option)
            self.in_model.parsamplfun(option);
            self.emis_model.parsamplfun(option);
            self.trans_model.parsamplfun(option);
        end
        function [out, sim_array, hmm2, modelstruct]=train(self,data,varargin)
            [out, sim_array, hmm2, modelstruct]=inference.train(self,data,varargin);
        end
        function [out modelstruct]=trainn(self,data,repi,varargin)
            [out modelstruct]=inference.trainn(self,data,repi,varargin);
        end
        function cleanpos(self)
           self.emis_model.cleanpos();
           self.in_model.cleanpos();
           self.trans_model.cleanpos();
        end
        function copy(self,hmm1)
            self.ndim=hmm1.ndim;
            self.nstates=hmm1.nstates;
            self.emis_model.copy(hmm1.emis_model);
            self.in_model.copy(hmm1.in_model);
            self.trans_model.copy(hmm1.trans_model);
        end
        function hmm2=copytrain(self,optcleanpos,nstates)
            if exist('nstates','var')
                newnstates=nstates;
            else
                newnstates=self.nstates;
            end
            pal=class(self.emis_model);
            e=eval([pal '(' num2str(self.ndim)  ')' ]);
            pal=class(self.in_model);
            i=eval([pal '(' num2str(newnstates)  ')' ]);
            pal=class(self.trans_model);
            t=eval([pal '(' num2str(newnstates)  ')' ]);
            hmm2=hmm(self.ndim,newnstates,e,i,t);
            
            %hmm2=hmm(self.ndim,newnstates);
            hmm2.copy(self);
            if optcleanpos==1
                hmm2.cleanpos();
            end
            hmm2.nstates=newnstates;
            hmm2.emis_model.nstates=newnstates;
            hmm2.trans_model.ndim=newnstates;
            hmm2.in_model.ndim=newnstates;
        end 
        function output=elibrefun(self,decodevar)
            self.emis_model.divklfun();
            output.divklobs=self.emis_model.divkl;
            self.in_model.divklfun();
            output.divklin=self.in_model.divkl;
            self.trans_model.divklfun();
            output.divkltrans=self.trans_model.divkl;
            output.loglik=decodevar.loglik;
            output.totaldivkl=-output.divkltrans-output.divklobs-output.divklin;
            output.fe=output.loglik+output.totaldivkl;
        end
        
        
    end
end

