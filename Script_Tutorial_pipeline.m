load Examples\Example5_EEG.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data file corresponds to EEG signal without artifacts (preprocessed) 
%and filtered in the interest frequency band (eg: 4- 30 Hz)
dmax=200;  %Maximun point in duration state
ncomp=10;   %PCA component
nstates=0;
forig=512; %Original Freq
fnew=128;  %New Freq
datahil=abs(hilbert(dataeeg64));  %Envolvent estimation
datares=resample(double(datahil),fnew,forig,2);
[coeff score]=pca(datares);
red=coeff(:,1:ncomp);
datain=datares*red;  %Data for model training
%%%%%MODEL CREATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e=emis_model.normal_normal_wishart(ncomp,nstates);
d=dur_model.lognormal_normal_gamma(1,nstates);
hsmm3=hsmm(ncomp,nstates,e,[],[],d);
hsmm3.priornoinf(datain,dmax);
hsmm3.emis_model.prior.prec_wishart{1}.scale=eye(ncomp)*mean(diag(inv(cov(datain))))/10/ncomp;
%Precision based in data
[outhsmm]=hsmm3.train(datain,'initoption','random','dmax',dmax,'maxitersim',1); %Model training
nstates2=outhsmm.nstates;
%%%%Topography recovering using GLM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y = zscore(datares); 
Xmat= zeros(length(outhsmm.stateseq),nstates2);
for j=1:nstates2
    Xmat(outhsmm.stateseq==j,j) = 1;
end
[beta] = util.spm.spm_robust_glm(Y, Xmat, 1,[]);
figure
for j=1:nstates2
    subplot(1,nstates2,j); 
    topoplot(beta(j,:), chanlocs, 'electrodes','off');     
end
