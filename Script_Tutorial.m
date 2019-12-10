%%%%%%%%%%%%% EXAMPLE 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state number unknown / default models;%%%%%%%%%%%%%%%%%%%%%%%%%%
load ./Examples/Example1_data2dim.mat  %Charge data, stateseq and data
ndim=2; %   Data dimension
nstates=0;   	%State number is unknown
hsmmmodel=hsmm(ndim,nstates);  %Create model (default option)
outhsmm=hsmmmodel.train(data); %Train model
%outhsmm have the simulation result
%hsmmmodel have the trained model
%%%%%%Most probably sequence%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq=hsmmmodel.viterbi(data);
%%%%%%Error estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq2=util.renameseq(stateseq,seq); %Rename of state using stateseq
error=sum(seq2~=stateseq)/length(data)*100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% EXAMPLE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state number known / change models;%%%%%%%%%%%%%%
load ./Examples/Example1_data2dim.mat  %Charge data, stateseq and data			
nstates=4;   	%State number is known
e=emis_model.normal_normal_gamma(ndim,nstates); %other Emision model
d=dur_model.lognormal_normal_gamma(1,nstates); %other Duration model
hsmmmodel=hsmm(ndim,nstates,e,[],[],d);  
outhsmm=hsmmmodel.train(data); %Train model
%%%%%%Most probably sequence%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq=hsmmmodel.viterbi(data);
%%%%%%Error estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq2=util.renameseq(stateseq,seq); %Rename of state using stateseq
error=sum(seq2~=stateseq)/length(data)*100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% EXAMPLE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state number known / change priors%%%%%%%%%%%%%%
load ./Examples/Example1_data2dim.mat  %Charge data, stateseq and data			
nstates=4;   	%State number is known
hsmmmodel=hsmm(ndim,nstates);  
hsmmmodel.priornoinf(data,100);  %Add non informative prior
hsmmmodel.dur_model.prior.mean_normal{1}.mean
hsmmmodel.dur_model.prior.mean_normal{1}.mean=40; %Change prior parameters
outhsmm=hsmmmodel.train(data); %Train model
%%%%%%Most probably sequence%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq=hsmmmodel.viterbi(data);
%%%%%%Error estimation%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seq2=util.renameseq(stateseq,seq); %Rename of state using stateseq
error=sum(seq2~=stateseq)/length(data)*100
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%RESULT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%See parameters posterior:
hsmmmodel.emis_model.posterior
hsmmmodel.trans_model.posterior
hsmmmodel.in_model.posterior
hsmmmodel.dur_model.posterior
%See paramters expectation:
hsmmmodel.emis_model.expectation 
hsmmmodel.trans_model.expectation  
hsmmmodel.in_model.expectation
hsmmmodel.dur_model.expectation
%See free energy curve
plot([outhsmm.fedetail.fe])
%See area of gamma
area(outhsmm.decodevar.gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%Change other options%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ./Examples/Example1_data2dim.mat  %Charge data, stateseq and data			
nstates=4;   	%State number is known
hsmmmodel=hsmm(ndim,nstates);
%Defaults Options:
%dmax=100;              %Maximun Duration
%initoption='random';   %Initialition (Random / Kmeans)
%maxcyc=100;            %Iteration Number in a run 
%maxitersim=1;          %Number of Run
%tol=0.01;              %out condition (% tolerance)
%minstates=2;           %State number searching begin in minstates
%maxstates=15;          %State number searching finish in maxtates 

%Example
outhsmm=hsmmmodel.train(data,'maxitersim',5,'dmax',200);











