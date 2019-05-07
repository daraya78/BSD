function ocupanciagraph(outhsmm,tit)
x0=100;
y0=50;
dx=200;
dy=200;
%f=figure;
%f.Position=[x0 y0  x0+dx y0+dy];
nstates=length(unique(outhsmm.stateseq));

[a1 b11]=hist(outhsmm.stateseq,1:nstates);
b=bar(b11,a1/length(outhsmm.stateseq)*100);
x1=xlabel('BS','FontSize',10);
y1=ylabel('BS Occupancy [%]','FontSize',10);
y1.FontAngle='italic';
x1.FontAngle='italic';
grid
title(tit)
%axis([0 nstates+1 0 100]);
end


