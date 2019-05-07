ndim=2;
nstates=4;
%Condicion 1 - Sujeto 1
e=emis_model.normal_normal_gamma(ndim,nstates);
hsmm3=hsmm(ndim,nstates,e); %Crea modelo
hsmm3.init(2); %Genera un sample de parametros para generar
data.cond(1).subj(1).block(1).data=hsmm3.gen(1000);
data.cond(1).subj(1).block(2).data=hsmm3.gen(800);
%Condicion 1 - Sujeto 2
hsmm3.init(2); %Genera un sample de parametros para generar
data.cond(1).subj(2).block(1).data=hsmm3.gen(1500);
data.cond(1).subj(2).block(2).data=hsmm3.gen(1500);
%Condicion 2 - Sujeto 1
hsmm3.init(2); %Genera un sample de parametros para generar
data.cond(2).subj(1).block(1).data=hsmm3.gen(1500);
%Condicion 2 - Sujeto 2
hsmm3.init(2); %Genera un sample de parametros para generar
data.cond(2).subj(2).block(1).data=hsmm3.gen(1000);
data.cond(2).subj(2).block(2).data=hsmm3.gen(1400);
data.cond(2).subj(2).block(3).data=hsmm3.gen(1400);










ncluster=10;
hmm3=hmm(ndim,ncluster); %Crea modelo
%Ingreso parametro para generar
hmm3.in_model.parsampl=ones(1,ncluster)/ncluster; %Condición inicial
hmm3.trans_model.parsampl=(ones(ncluster)-eye(ncluster))/(ncluster-1);%Transición
%Emisión
for j=1:ncluster
    hmm3.emis_model.parsampl.mean{j}=rand(1,ndim)*20; %Cambiando el 20 por uno menor hay mas traslape
    hmm3.emis_model.parsampl.prec{j}=ones(1,ndim);
end

%Generar
[data seq]=hmm3.gen(1000);

%Graficar
figure
for j=1:ncluster
    plot(data(seq==j,1),data(seq==j,2),'.')
    hold on
end


























