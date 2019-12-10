opt.train='VB';
opt.maxcyc=100;%100;
opt.tol=0.1;
opt.dmax=100;
%Kmeans
opt.nrep=10;
opt.maxiter=50000;
opt.initoption='random';
%%%%%%%%%%%%%%%%%%
opt.maxitersim=1;
opt.minstates=2;
opt.maxstates=15;
opt.seq='maxgamma';
opt.verbose='yes';
opt.prepro='nothing';
opt.decalpha=1;
opt.decbeta=1;
opt.deceta=0;
opt.prior='databased';
opt.shareemis = [1 0 1];
opt.sharetrans = [1 0 1];
opt.sharein = [0 1 1];
opt.sharedur = [1 0 1];
opt.parallel=0;

