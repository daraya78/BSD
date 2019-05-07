ndim=2;
nstates=3;

%%%CASO TODO DIFIRENTE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Condicion 1 - Sujeto 1
e=emis_model.normal_normal_wishart(ndim,nstates);
hsmm3=hsmm(ndim,nstates,e); %Crea modelo
hsmm3.init(2); %Genera un sample de parametros para generar
[data.cond(1).subj(1).block(1).data seq{1,1,1}]=hsmm3.gen(1000);
[data.cond(1).subj(1).block(2).data seq{1,1,2}]=hsmm3.gen(800);
%Condicion 1 - Sujeto 2
hsmm3.init(2); %Genera un sample de parametros para generar
[data.cond(1).subj(2).block(1).data seq{1,2,1}]=hsmm3.gen(1500);
[data.cond(1).subj(2).block(2).data seq{1,2,2}]=hsmm3.gen(1500);
%Condicion 2 - Sujeto 1
hsmm3.init(2); %Genera un sample de parametros para generar
[data.cond(2).subj(1).block(1).data seq{2,1,1}]=hsmm3.gen(1500);
%Condicion 2 - Sujeto 2
hsmm3.init(2); %Genera un sample de parametros para generar
[data.cond(2).subj(2).block(1).data seq{2,2,1}]=hsmm3.gen(1000);
[data.cond(2).subj(2).block(2).data seq{2,2,2}]=hsmm3.gen(1400);
[data.cond(2).subj(2).block(3).data seq{2,2,3}]=hsmm3.gen(1400);

%%%CASO COMPARTE MODELO DE EMISIÓN PERO NO DE DURACIÓN%%%%%%%%%%%%%%%%%%%%

e1=emis_model.normal_normal_wishart(ndim,nstates);
e2=emis_model.normal_normal_wishart(ndim,nstates);
e1.parsamplfun(2);
e2.copy(e1);
hsmm3=hsmm(ndim,nstates); %Crea modelo
%Condicion 1 - Sujeto 1
hsmm3.init(2);
hsmm3.emis_model=e2;
[data.cond(1).subj(1).block(1).data seq{1,1,1}]=hsmm3.gen(2000);
[data.cond(1).subj(1).block(2).data seq{1,1,2}]=hsmm3.gen(2800);
%Condicion 1 - Sujeto 2
hsmm3.init(2); %Genera un sample de parametros para generar
e2.copy(e1);
hsmm3.emis_model=e2;
[data.cond(1).subj(2).block(1).data seq{1,2,1}]=hsmm3.gen(2500);
[data.cond(1).subj(2).block(2).data seq{1,2,2}]=hsmm3.gen(2500);
%Condicion 2 - Sujeto 1
hsmm3.init(2); %Genera un sample de parametros para generar
e2.copy(e1);
hsmm3.emis_model=e2;
[data.cond(2).subj(1).block(1).data seq{2,1,1}]=hsmm3.gen(2500);
%Condicion 2 - Sujeto 2
hsmm3.init(2); %Genera un sample de parametros para generar
e2.copy(e1);
hsmm3.emis_model=e2;
[data.cond(2).subj(2).block(1).data seq{2,2,1}]=hsmm3.gen(2000);
[data.cond(2).subj(2).block(2).data seq{2,2,2}]=hsmm3.gen(2400);
[data.cond(2).subj(2).block(3).data seq{2,2,3}]=hsmm3.gen(2400);













%Graficar
figure
conta=1;
for kcond=1:2
    for ksubj=1:2
        nblock=size(data.cond(kcond).subj(ksubj).block,2);
        subplot(2,2,conta)
        for kblock=1:nblock
            data2=data.cond(kcond).subj(ksubj).block(kblock).data;
            for j=1:nstates
                plot(data2(seq{kcond,ksubj,kblock}==j,1),data2(seq{kcond,ksubj,kblock}==j,2),'x')
                hold on
            end    
        end
        title(['Cond ' num2str(kcond) ' Subj: ' num2str(ksubj)])
        conta=conta+1;
    end
end























