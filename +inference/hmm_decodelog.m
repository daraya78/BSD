        function [decodevar] = hmm_decodelog(self,X,opt_varargin)
            n=0;            
            if isempty(opt_varargin)
                program_option;
                n=size(opt_varargin,2);
            elseif isstruct(opt_varargin{1})
                opt=opt_varargin{1};
            else
                program_option;
                n=size(opt_varargin,2);
            end
            for j=1:2:n
                opt=setfield(opt,opt_varargin{j},opt_varargin{j+1});
            end
            
            [B2 B]  = self.emis_model.prob(opt.train,X);
            [A2 A]  = self.trans_model.expect(opt.train);
            [in2 in] = self.in_model.expect(opt.train);
            T=size(B,1);
            alpha(1,:)=in+B(1,:);
            scale=-inf(T,1);
            scale(1)=util.logsumexp(alpha(1,:),2);
            alpha(1,:)=alpha(1,:)-scale(1);
            for i=2:T
                alpha(i,:)=(util.logmultmat(alpha(i-1,:),A))+B(i,:);
                scale(i)=util.logsumexp(alpha(i,:),2);
                alpha(i,:)=alpha(i,:)-scale(i);
            end;
            beta(T,:)=zeros(1,self.nstates)-scale(T);
            for i=T-1:-1:1
                beta(i,:)=util.logmultmat((beta(i+1,:)+B(i+1,:)),A')-scale(i);
            end;
            gamma=(alpha+beta);
            gamma=bsxfun(@minus,gamma,util.logsumexp(gamma,2));
            pSeq = sum(scale);
            xi=-inf(self.nstates,self.nstates,T-1);
            for i=1:T-1
                t=A+util.logmultmat(alpha(i,:)',(beta(i+1,:)+B(i+1,:)))-scale(i+1);
                xi(:,:,i)=t-util.logsumexp(t(:));
            end;
            decodevar.alpha=exp(alpha);
            decodevar.beta=exp(beta);
            decodevar.Emi=B2;
            decodevar.scale=exp(scale);
            decodevar.gamma=exp(gamma);
            decodevar.xi=exp(xi);
            decodevar.loglik=pSeq;
        end