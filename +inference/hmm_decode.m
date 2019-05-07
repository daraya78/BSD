        function [decodevar] = hmm_decode(self,X,opt_varargin)
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
            B  = self.emis_model.prob(opt.train,X);
            A  = self.trans_model.expect(opt.train);
            in = self.in_model.expect(opt.train);
            T=size(B,1);
            alpha(1,:)=in.*B(1,:);
            scale=zeros(T,1);
            scale(1)=sum(alpha(1,:));
            alpha(1,:)=alpha(1,:)/scale(1);
            for i=2:T
                alpha(i,:)=(alpha(i-1,:)*A).*B(i,:);
                scale(i)=sum(alpha(i,:));
                alpha(i,:)=alpha(i,:)/scale(i);
            end;
            beta(T,:)=ones(1,self.nstates)/scale(T);
            for i=T-1:-1:1
                beta(i,:)=(beta(i+1,:).*B(i+1,:))*(A')/scale(i);
            end;
            gamma=(alpha.*beta);
            gamma=bsxfun(@rdivide,gamma,sum(gamma,2));
            pSeq = sum(log(scale));
            xi=zeros(self.nstates,self.nstates,T-1);
            for i=1:T-1
                t=A.*( alpha(i,:)' * (beta(i+1,:).*B(i+1,:)))/scale(i+1);
                xi(:,:,i)=t./sum(t(:));
            end;
            decodevar.alpha=alpha;
            decodevar.beta=beta;
            decodevar.Emi=B;
            decodevar.scale=scale;
            decodevar.gamma=gamma;
            decodevar.xi=xi;
            decodevar.loglik=pSeq;
        end