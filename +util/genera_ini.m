function [ini]=genera_ini (p,ndata,dmax,nite,nstatesmax)
ruta0=cd;
rutapro='P:\BrainDynamicsLab\Software\HSMM\Software HSMM\hsmm31Mar2020';
cd(rutapro)
for kstate=2:nstatesmax
    for j=1:nite
        ini{kstate}(j,:)=util.unirndseq(ndata-p,kstate,dmax);
    end
end
cd(ruta0)
end







