opt.train='VB';
opt.maxcyc=100;%100;
opt.tol=0.01;
opt.dmax=400;
%Kmeans
opt.nrep=10;
opt.maxiter=50000;
opt.initoption='kmeans';
%%%%%%%%%%%%%%%%%%
opt.maxitersim=1;
opt.minstates=1;
opt.maxstates=15;
opt.seq='maxgamma';
opt.verbose='yes';
opt.prepro='nothing';
opt.decalpha=1;
opt.decbeta=1;
opt.deceta=0;
opt.prioremis='default';
opt.shareemis = 1;
opt.sharetrans = 0;
opt.sharecond = 0;
opt.sharedur = 0;
opt.shareemis = [1 0 1];
opt.sharetrans = [1 0 1];
opt.sharein = [0 1 1];
opt.sharedur = [1 0 1];

