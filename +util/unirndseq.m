function re=unirndseq(ndata,nstates,dmax);
    re=[];
    statek1=randi(nstates);
    re=repmat(statek1,1,randi(dmax));
    while (size(re,2)<ndata)
        statek2=randi(nstates);
        rang=ceil([dmax/5-dmax/8 dmax/5+dmax/8]);
        if statek2~=statek1
            %re=[re repmat(statek2,1,randi(dmax))];
            aux=randi(rang);
            re=[re repmat(statek2,1,aux)];
            statek1=statek2;
        end
    end
    re=re(1:ndata);
end
