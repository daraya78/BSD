        function [decodevar] = hsmmresidual_decodelog(self,X,opt_varargin)
            
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
            [B2 B]  = self.emis_model.prob(opt.train,X);
            [A2 A]  = self.trans_model.expect(opt.train);
            [in2 in] = self.in_model.expect(opt.train);
            [P2 P]  = self.dur_model.prob(opt.train,(1:1:dmax)');
            T=size(B,1);
            B=B';
            P=P';
            m=self.nstates;
            %ALPHA=[in', zeros(m, dmax-1)];
            ALPHA=repmat(in',1,dmax);
            
            for t=1:T
                x=repmat((util.logmultmat(A',ALPHA(:,1))+B(:,t)),1,dmax)+P;
                w=ALPHA(:,2:dmax)+repmat(B(:,t),1,dmax-1);
                ALPHA=util.logsumexpmat([w,-inf(m,1)],x);
                c(t)=-util.logsumexp(ALPHA(:));
                ALPHA=ALPHA+c(t);
                ALPHAx(:,t)=ALPHA(:,1);
                if opt.decalpha==1
                    ALPHATOTAL(:,:,t)=ALPHA;
                end
            end
            %BETA=[ones(m,1),zeros(m,dmax-1)];
            %BETA=repmat(zeros(m,1),1,dmax);  
            BETA=[zeros(m,1),-inf(m,dmax-1)];
            
            
            GAMMA=util.logsumexp(ALPHA,2);     
            GAMMATOTAL(T,:)=GAMMA;
            if opt.decbeta==1
                BETATOTAL(:,:,T)=BETA;
            end
            durcount=[];
            for t=(T-1):-1:1
                y=B(:,t+1)+c(t+1);
                z=y+(util.logsumexp((P+BETA),2));
                zdur=repmat(y,1,dmax)+P+BETA;
                BETA(:,2:dmax)=repmat(y,1,dmax-1)+BETA(:,1:dmax-1);
                BETA(:,1)=util.logmultmat(A,z);
                if opt.decbeta==1
                    BETATOTAL(:,:,t)=BETA;
                end
                aux1=repmat(ALPHAx(:,t),[1 dmax]);
                aux2=repelem(reshape(aux1,[m*dmax 1]),m);
                aux3=reshape(repelem(zdur,1,m),[m*m*dmax 1]);
                aux4=sum([aux2 aux3 repmat(reshape(A',[m*m 1]),[dmax 1])],2);
                aux5=reshape(aux4,[m m dmax]);
                aux6=permute(aux5,[2 1 3]);
                aux7=squeeze(util.logsumexp(aux6,1));
                aux8=permute(aux7,[2 1 3]);
                if opt.deceta==1
                    eta(:,:,t)=aux8;
                end
                durcount=squeeze(util.logsumexp(cat(3,aux8,durcount),3));
                XI=(util.logmultmat(ALPHAx(:,t),z'))+A;
                XITOTAL(:,:,t)=XI;
                XIaux1=util.logsumexpmat(GAMMA,util.logsumexp(XI,2));
                XIaux2=util.logsumexp(XI',2);
                %Temporalmente se modifica
                %GAMMA=util.logsubexpmat(XIaux1,XIaux2);
                %GAMMATOTAL(t,:)=real(GAMMA);
            end
            
            %Cambio provisorio
            %gamma=GAMMATOTAL;
            %nor=util.logsumexp(gamma(10,:),2);
            %gamma=gamma-nor;
            aux=ALPHATOTAL+BETATOTAL;
            gamma=squeeze(sum(exp(aux),2))';
            nor=sum(gamma(10,:));
            gamma=gamma/nor;
            nor=log(nor);
            gamma=log(gamma);
            
            
            xi=XITOTAL-nor;
            durcount=durcount-nor;
            if opt.deceta==1
                eta=eta-nor;
            end
            loglik=sum(-c); 
            if isnan(loglik)
                'loglik indefinida'
            end
            %decodevar.Emi=exp(B2); 
            decodevar.scale=exp(c);
            decodevar.gamma=exp(gamma);
            decodevar.xi=exp(xi);
            decodevar.loglik=loglik;
            decodevar.durcount=exp(durcount);
            decodevar.alpha=[];
            decodevar.beta=[];
            decodevar.eta=[];
            if opt.decalpha==1
                %decodevar.alpha=exp(ALPHATOTAL);
            end
            if opt.decbeta==1
                %decodevar.beta=exp(BETATOTAL);
            end
            if opt.deceta==1
                %decodevar.eta=exp(eta);
            end
        end
        

