ncomp=10;
%nstates=4;
%kprec=1000000; % Incluso con este número tan alto no converge a un número
%pequeño de estados, se propone utilizar directamente nstates=4

%rutadata='P:\BrainDynamicsLab\Software\HSMM\Tutorial Nelson';
rutadata='P:\BrainDynamicsLab\Software\HSMM\Trabajo TMS';

rutaprog='P:\BrainDynamicsLab\Software\HSMM\hsmm_v19Mar2019';
forig=512;  %Frecuencia en la cual estan los datos
fnew=128;
cd(rutadata)
%load('cond3subj1block1.mat');
load('EEGsuj1.mat');
datahil1=abs(hilbert(EEG.data'));
datares1=resample(double(datahil1),fnew,forig,2);
datares2 = bsxfun(@minus, datares1, mean(datares1));
[coeff score]=pca(datares2);
red=coeff(:,1:ncomp);
data1=datares2*red;


[datablock1 tramo1 tramoelim1 type1]=segmenta(data1,EEG.event,forig,fnew);
%This function segment signal in base to the marks and create the
%datastruct 

cd(rutaprog)
e=emis_model.normal_normal_wishart(ncomp,nstates);
d=dur_model.lognormal_normal_gamma(1,nstates);
hsmm3=hsmm(ncomp,nstates,e,[],[],d);
%hsmm3.emis_model.prior.prec_wishart{1}.scale=diag(1./(((max(data1)-min(data1))/2).^2))/ncomp/kprec;

%Flag that determine that features is share (column are: cond / subj / block
shareemis =  [1 1 1];
sharetrans = [0 0 1];
sharein =    [0 0 1];
sharedur =   [0 0 1];

[outhsmm, sim_array, hsmmarray modelstruct]=hsmm3.train(datablock1,'initoption', ...
    'random','maxitersim',1,'shareemis',shareemis,'sharetrans',sharetrans,'sharein', ...
    sharein,'sharedur',sharedur);

%Estructura outhsmm contiene los resultados, en este caso incluye la seq de
%estados
%Estructura modelstruct contiene todos los modelos asociados a las
%condiciones, sujetos y bloques.


%%%%RECUPERACION DE TOPOGRAFIAS USANDO PCA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SE GRAFICA LAS TOPOGRAFIAS DEL PRIMER MODELO (TODAS SON IGUALES (SE
%COMPARTE MODELO DE EMISION)
f=figure;
f.Position=[1376 227 985 254];
nstates2=outhsmm.nstates.cond(1).subj(1).block(1).nstates;
for kstate=1:nstates2
    m1=modelstruct.matrixmodel(1).emis_model.expectation.mean{kstate};
    topoaux=m1*red';
    subplot(1,nstates2,kstate); 
    topoplot(topoaux, EEG.chanlocs, 'electrodes','off'); 
    title(['State ' num2str(kstate)])
end


%%%%RECUPERACION DE TOPOGRAFIAS USANDO GLM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(rutadata)
[datablock2 tramo1 tramoelim1 type1]=segmenta(datares2,EEG.event,forig,fnew);
nmax=min(20000,size(datablock2.cond(1).subj(1).block(1).data,1));
upoints=1:nmax;
Y = zscore(datablock2.cond(1).subj(1).block(1).data(upoints,:)); 
seq1=outhsmm.stateseq.cond(1).subj(1).block(1).stateseq;
seqest = seq1(upoints);
Xmat= zeros(length(seqest),nstates2);
for j=1:nstates2
    Xmat(seqest==j,j) = 1;
end
[beta, sigma, E, V ] = mvregress(Xmat,Y,'algorithm','cwls');
f=figure;
f.Position=[1399 256 966 208];
for kstate=1:nstates2
    subplot(1,nstates2,kstate); 
    topoplot(beta(kstate,:), EEG.chanlocs, 'electrodes','off');
    title(['State ' num2str(kstate)])
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Grafico de Ocupancia
f=figure;
f.Position=[1407 193 889 225];
for kcond=1:size(outhsmm.nstates.cond,2)
    subplot(1,size(outhsmm.nstates.cond,2),kcond)
    nstates2=outhsmm.nstates.cond(kcond).subj(1).block(1).nstates;
    stateseq=outhsmm.stateseq.cond(kcond).subj(1).block(1).stateseq;
    [a b]=hist(stateseq,1:nstates2);
    b=bar(b,a/length(stateseq)*100);
    axis([0 nstates+1 0 50])
    grid
    title(['Condition ' num2str(kcond)])
    
end

%Grafico de matriz de transicion en cada condicion
f=figure;
f.Position=[66 281 1151 236];
for kcond=1:size(outhsmm.nstates.cond,2)
    subplot(2,ceil(size(outhsmm.nstates.cond,2)/2),kcond)
    matriz=modelstruct.matrixmodel(kcond,1,1).trans_model.expectation;
    imagesc(matriz)
    colorbar
    title(['Condition ' num2str(kcond)])
end

%Grafico de pdf de duración


f=figure;
f.Position=[24 303 1241 288];
t=0:0.1:100;
for kcond=1:size(outhsmm.nstates.cond,2)
    subplot(2,ceil(size(outhsmm.nstates.cond,2)/2),kcond)
    nstates2=outhsmm.nstates.cond(kcond).subj(1).block(1).nstates;
    for kstate=1:nstates2
        m2=modelstruct.matrixmodel(kcond,1,1).dur_model.expectation.mean{kstate};
        p2=modelstruct.matrixmodel(kcond,1,1).dur_model.expectation.prec{kstate};
        plot(t,lognpdf(t,m2,inv(p2)))
        hold on
        le{kstate}=['State ' num2str(kstate)];
    end
    grid
    legend(le)
    title(['Condition ' num2str(kcond)])
    axis([0 50 0 0.2])
end










