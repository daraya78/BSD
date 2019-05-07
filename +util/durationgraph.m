function durationgraph(hsmm3,srate,tit)

x0=100;
y0=50;
dx=400;
dy=300;

xtime=0.1:0.1:200;
dis1=hsmm3.dur_model.prob('VB',xtime');

%f=figure;
%f.Position=[x0 y0  x0+dx y0+dy];
ttime=200*1/srate;
lxtime=length(xtime);
t=0:ttime/lxtime:ttime-ttime/lxtime;
p1=plot(t,dis1);
grid
ti=title(tit,'FontSize',12);
%ti.Position=[670.0009 0.0325 0];
ylabel('q(d)','FontSize',10)
xlabel('d [s]','FontSize',10)
s1.Position=[0.1700 0.6378 0.7750 0.2872];
s1.YLabel.FontAngle='italic';

%te1=text(500,0.018,'Log-Normal','FontSize',13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:hsmm3.nstates
   le{j}=['Estado: ' num2str(j)]; 
end
legend(le)

end