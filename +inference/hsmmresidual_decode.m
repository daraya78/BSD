        function [decodevar] = hsmmresidual_decode(self,X,opt_varargin)

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
            ma=[];
            for j=1:2:n
                opt=setfield(opt,opt_varargin{j},opt_varargin{j+1});
            end
        
            dmax=opt.dmax;
            B  = self.emis_model.prob(opt.train,X);
            A  = self.trans_model.expect(opt.train);
            in = self.in_model.expect(opt.train);
            P  = self.dur_model.prob(opt.train,(1:1:dmax)');
            T=size(B,1);
            B=B';
            P=P';
            m=self.nstates;
            %ALPHA=[in', zeros(m, dmax-1)];
            ALPHA=repmat(in',1,dmax);
            
            for t=1:T
                x=repmat(((A'*ALPHA(:,1)).*B(:,t)),1,dmax).*P;
                w=ALPHA(:,2:dmax).*repmat(B(:,t),1,dmax-1);
                ALPHA=[w,zeros(m,1)]+x;
                c(t)=1/sum(ALPHA(:)); 
                ALPHA=ALPHA.*c(t); 
                ALPHAx(:,t)=ALPHA(:,1);
                ALPHATOTAL(:,:,t)=ALPHA;
            end
            
            %BETA=[ones(m,1),zeros(m,dmax-1)];
            BETA=repmat(ones(m,1),1,dmax);
            
            BETATOTAL(:,:,T)=BETA;
            GAMMA=sum(ALPHA,2);
            [u,S_est(T)]=max(GAMMA);
            GAMMATOTAL(T,:)=GAMMA;
            for t=(T-1):-1:1
                y=B(:,t+1).*c(t+1);
                z=y.*(sum((P.*BETA),2));
                zdur=repmat(y,1,dmax).*P.*BETA;
                BETA(:,2:dmax)=repmat(y,1,dmax-1).*BETA(:,1:dmax-1);
                BETA(:,1)=A*z;
                BETATOTAL(:,:,t)=BETA;
                for d=1:dmax
                    XIdur(:,:,d,t)=(ALPHAx(:,t)*zdur(:,d)'.*A);
                end
                XI=(ALPHAx(:,t)*z').*A;
                XITOTAL(:,:,t)=XI;
                GAMMA=GAMMA+sum(XI,2)-sum(XI',2);
                GAMMATOTAL(t,:)=GAMMA;
                [u,S_est(t)]=max(GAMMA);
            end

            aux=ALPHATOTAL.*BETATOTAL;
            aux2=sum(aux,2);
            for t=1:T
                gamma(t,:)=aux2(:,1,t)';
            end    
            nor=sum(gamma(10,:));
            gamma=gamma/nor;
            XITOTAL=XITOTAL/nor;
            clear aux
            aux3=sum(XIdur,1);
            for t=1:T-1
                aux(:,:)=aux3(1,:,:,t);
                eta(:,:,t)=aux';
            end
            xi=XITOTAL;
            eta=eta/nor;
            clear aux;
            clear aux2;
            aux=sum(sum(XIdur,1),4);
            durcount(:,:)=aux(1,:,:)/nor;
            
            loglik=sum(log(1./c));
            if isnan(loglik)
                'loglik indefinida'
            end
            
            
           % decodevar.alpha2=ALPHA2;
            decodevar.Emi=B; 
            decodevar.alpha=ALPHATOTAL;
            decodevar.beta=BETATOTAL;
            decodevar.scale=c;
            decodevar.gamma=gamma;
            decodevar.eta=eta;
            decodevar.xi=xi;
            decodevar.loglik=loglik;

        end
        

