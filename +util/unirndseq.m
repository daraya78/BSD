function re=unirndseq(ndata,nstates,dmax);
    re=[];
    rang=ceil([dmax/5-dmax/8 dmax/5+dmax/8]);
    statek1=randi(nstates);
    re=repmat(statek1,1,randi(rang));
    while (size(re,2)<ndata)
        statek2=randi(nstates);
        
        if statek2~=statek1
            %re=[re repmat(statek2,1,randi(dmax))];
            re=[re repmat(statek2,1,randi(rang))];
            statek1=statek2;
        end
    end
    re=re(1:ndata);
end
