function multi2=singletomulti(single,multi)


hsmm5=hsmm();
multi2=multimodel(hsmm5);
[ny nx]=size(multi.matrixmodel);
multi2.n=multi.n;
multi2.iall=multi.iall;
multi2.iemis=multi.iemis;
multi2.itrans=multi.itrans;
multi2.iin=multi.iin;
multi2.idur=multi.idur;

for ky=1:ny
    for kx=1:nx
        d=dur_model.lognormal_normal_gamma();
        multi2.matrixmodel(ky,kx)=hsmm([],[],[],[],[],d);
        multi2.matrixmodel(ky,kx).copy(single)
    end
end

end




