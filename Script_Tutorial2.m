ruta1='P:\BrainDynamicsLab\Software\HSMM\Sim_paper\Realdata_13Nov\';
ruta2='P:\BrainDynamicsLab\Software\HSMM\hsmm_v19Mar2019\';
file='s7b1 clean broadband ENV.mat';
file2='s7b1 clean broadband.mat';
load([ruta1 file]);
load([ruta1 file2]);
data=ADATA';
datalog=log(data-min(data)+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data file corresponds to EEG signal without artifacts (preprocessed) 
%and filtered in the interest frequency band (eg: 4- 30 Hz)
ncomp=5;
nstates=0;
kprec=1000;
forig=512; %Original Freq
fnew=128;  %New Freq
datahil1=abs(hilbert(datalog));
datares1=resample(double(datahil1),fnew,forig,2);
datares2 = bsxfun(@minus, datares1, mean(datares1));
[coeff score]=pca(datares2);
red=coeff(:,1:ncomp);
data1=datares2*red;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e=emis_model.normal_normal_wishart(ncomp,nstates);
d=dur_model.lognormal_normal_gamma(1,nstates);
hsmm3=hsmm(ncomp,nstates,e,[],[],d);
%hsmm3.emis_model.prior.prec_wishart{1}.scale=diag(diag(inv(cov(data1))))/ncomp/kprec;
[outhsmm]=hsmm3.train(data1,'initoption','random','minstates',5,'maxstates',10,'maxitersim',1);
%seq=hsmmm3.viterbi(data1);
nstates2=outhsmm.nstates;
%%%%RECUPERACION DE TOPOGRAFIAS USANDO GLM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
upoints=1:20000;
Y = zscore(data2(upoints,:)); 
seqest = seq(upoints);
Xmat= zeros(length(seqest),nstates2);
for j=1:nstates2
    Xmat(seqest==j,j) = 1;
end
[beta, sigma, E, V ] = mvregress(Xmat,Y,'algorithm','cwls');
figure
for j=1:nstates2
    subplot(1,nstates2,j); 
    topoplot(beta(j,:), EEG.chanlocs, 'electrodes','off');     
end
%%%%RECUPERACION DE TOPOGRAFIAS USANDO PCA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:nstates2
    topoaux=hsmmmodel.emis_model.posterior.mean_normal{1}.mean*red';
    subplot(1,nstates2,j); 
     topoplot(topoaux, EEG.chanlocs, 'electrodes','off'); 
end
