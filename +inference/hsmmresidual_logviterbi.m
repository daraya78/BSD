function [stateseq delta sa] = hsmmresidual_logviterbi(self,X,opt_varargin)
    nstates=self.nstates;
    ndata=self.emis_model.ndatafun(X);
    
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
    [B2 B3]  = self.emis_model.prob(opt.train,X);
    [A2 A]  = self.trans_model.expect(opt.train);
    [in2 ini] = self.in_model.expect(opt.train);
    [P2 pdur]  = self.dur_model.prob(opt.train,(1:1:dmax)');

       
%    for j=1:self.nstates
%        EMI(:,j) = self.obs_model.prob(2,self.data,j);
%        pdur(:,j) = self.dur_model{j}.durprob(2,1:1:dmax);
%    end
%    ini=self.in_model.inprob(2);    
    for j=1:self.nstates
        B(j,:,:) = util.pdfdur(B3(:,j),dmax); %create u(state,t,dur)
    end
       
    delta=-Inf(nstates,ndata,dmax);
    deltaaux=-Inf(nstates,ndata,dmax);
    
    delta(:,1,1)=ini;

    
    %FORDWARD ALGORITHM
    for t=1:ndata+dmax-1
        if mod(t,1000)==0
            t
        end
        for j=1:nstates
            for d=1:dmax
                deltamax=-inf;
                imax=0;
                hmax=0;
                %for i=1:nstates
                %    for h=1:dmax
                        if (t-d)==0
                            %deltaaux(j,t,d)=ini(1,j)+A(i,j)+pdur(d,j)+B(j,t,d);
                            deltaaux=ini(1,j)+A(:,j)+pdur(d,j)+B(j,t,d);
                        elseif (t-d)>=1 & t<=ndata
                            %deltaaux(j,t,d)=delta(i,t-d,h)+pdur(d,j)+A(i,j)+B(j,t,d);

                            try
                            deltaaux=squeeze(delta(:,t-d,:))+pdur(d,j)+repmat(A(:,j),1,dmax)+B(j,t,d);
                            catch
                                keyboard
                            end
                            
                            
                        else
                            deltaaux=-inf;
                            posx=0;
                            posy=0;
                        end
                        
                        if (t-d)>=0 & t<=ndata
                            %if deltaaux(j,t,d)>deltamax
                            %    delta(j,t,d)=deltaaux(j,t,d);
                            %    deltamax=delta(j,t,d);
                            %    imax=i;
                            %    hmax=h; 
                            %end
                            
                            [m py]=max(deltaaux);
                            [m posx]=max(m);
                            posy=py(posx);
                            delta(j,t,d)=deltaaux(posy,posx);
                            if isinf(deltaaux(posy,posx))
                                posy=0;
                                posx=0;
                            end
                            
                            
                            
                        end                           
                %    end
                %end
                %sa(t,j,d,:)=[t-d imax hmax];
                sa(t,j,d,:)=[t-d posy posx];
                
            end
        end
    end
  %THE TRACE BACK ALGORITHM
  finalstates=-inf;
  for j=1:nstates
      for d=1:dmax
          if delta(j,ndata,d) > finalstates
              finalstates=delta(j,ndata,d);
              jm=j;
              dm=d;
          end
      end
  end
  stateseq=ones(1,dm)*jm;
  count=ndata-dm;
  while (count & dm) > 0
      dmaux=sa(count+dm,jm,dm,3);
      jmaux=sa(count+dm,jm,dm,2);
      stateseq=[ones(1,dmaux)*jmaux stateseq];
      jm=jmaux;
      dm=dmaux;
      count=count-dm;
  end
  
  
end