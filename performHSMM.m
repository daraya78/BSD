%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PERFORM HSMM Multi subject
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% SET HSMM OPTIONS
% Case 2:                                    		
%                                    Condition	Subject	Block
% To Share Emission Model:  	   	 SI	        SI		SI     
% To share transition model: 		 NO			NO		SI
% To share duration model:           NO			NO		SI
% To share initial condition model:  NO			NO		SI

ndim = 20; % number of retained components for pca
nstates = 4; % define number of states to model

% EMISSION MODEL: NORMAL DIST WITH FULL COV MATRIX
% ---------------------------
e = emis_model.normal_normal_wishart(ndim,nstates);
hsmm3 = hsmm(ndim,nstates);
hsmm3.emis_model = e;

% define shared models for multi subject HSMM analysis
se = [1 1 1]; % Emission model SHARED between cond,subj and block
st = [0 1 1]; % transition model NOT SHARED between cond,subj and SHARED between block
si = [0 0 1]; % initial state model NOT SHARED between cond,subj and SHARED between block
sd = [0 0 1]; % duration model NOT SHARED between cond,subj and SHARED between block

% Emission model Prior:
% ---------------------------
% Mean: By default it will have mean prior: % Mean=0; % Prec=0.1 (diag matrix)
hsmm3.emis_model.prior.mean_normal{1}.mean = zeros(1,ndim);
hsmm3.emis_model.prior.mean_normal{1}.prec = eye(ndim)*0.1;
% i.e., it should be check that data after PCA have values into this range

% Precision: By default it will have precision prior:% Shape=0.001;% Scale=1000
% This is a very common non-informative prior 
hsmm3.emis_model.prior.prec_wishart{1}.scale = eye(ndim)*0.01;
hsmm3.emis_model.prior.prec_wishart{1}.degree = 20;

% DURATION MODEL: LOG NORMAL
% ---------------------------
hsmm3.dur_model = dur_model.lognormal_normal_gamma(1,nstates); % corrected
%hsmm3.dur_model = dur_model.normal_normal_gamma(1,nstates); % corrected
% Duration model Prior: 
% Mean:
hsmm3.dur_model.prior.mean_normal{1}.mean = 2;
hsmm3.dur_model.prior.mean_normal{1}.prec = 0.1;%0.2;
% Precision:
hsmm3.dur_model.prior.prec_gamma{1}.shape = 0.001;
hsmm3.dur_model.prior.prec_gamma{1}.scale = 1000;






%% TRAINING
dmax = 100; %20; %200;   % Maximum duration in points
nrep = 10;	  % Number of times to initialize (e.g. Kmeans)
initoption = 'random';	% Initialization algorithm
maxitersim = 10; %5% N times to repeat the algorithm (the final HSMM will be the one with the max free energy)
% minstates = 2; % Minimum number of states
% maxstates = 6; % Maximum number of states





tic;
% Use the training command indicating the options:
[out3, aux, aux2, matrix] = hsmm3.train(data1,'dmax',dmax,'nrep',nrep,...
    'initoption',initoption,'maxitersim',maxitersim,'shareemis',...
    se,'sharetrans',st,'sharein',si,'sharedur',sd);
% outdatahsmm20 = hsmm20.train(data,'dmax',dmax,'nrep',nrep,'initoption',initoption,'maxitersim',maxitersim);
%,'minstates',minstates,'maxstates',maxstates);
ttrain = toc;
disp('Training finished');



data=data1.cond(1).subj(1).block(1).data;
[out3, aux, aux2] = hsmm3.train(data,'dmax',dmax,'nrep',nrep,...
    'initoption',initoption,'maxitersim',maxitersim,'shareemis',...
    se,'sharetrans',st,'sharein',si,'sharedur',sd,'seq','maxgamma');

