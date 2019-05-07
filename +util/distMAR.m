function re=distMAR(a,b);
    
    if strcmp(class(m1),'varm')
        m1=a;
        m2=b;
        for k=1:30
            data=simulate(m1,1000);
            dim=size(data,2);
            resid1=infer(m1,data);
            resid2=infer(m2,data);
            nor=log(1/((2*pi)^(dim/2)*det(eye(2))^(1/2) ));
            loglik1(k)=sum(nor-sum((1/2*(resid1-zeros(1,2))*inv(eye(2)).*(resid1-zeros(1,2))),2));
            loglik2(k)=sum(nor-sum((1/2*(resid2-zeros(1,2))*inv(eye(2)).*(resid2-zeros(1,2))),2));
        end
        re=(mean(loglik1)-mean(loglik2))/1000;
    else %It is Coef
        for k=1:30
            %ret=sqrt(size(a,2)/ndim);
            %m1=varm(2,2);
            %m2=varm(2,2);
            %m1.AR{1}=
            %PENDIENTE
            data=simulate(m1,1000);
            dim=size(data,2);
            resid1=infer(m1,data);
            resid2=infer(m2,data);
            nor=log(1/((2*pi)^(dim/2)*det(eye(2))^(1/2) ));
            loglik1(k)=sum(nor-sum((1/2*(resid1-zeros(1,2))*inv(eye(2)).*(resid1-zeros(1,2))),2));
            loglik2(k)=sum(nor-sum((1/2*(resid2-zeros(1,2))*inv(eye(2)).*(resid2-zeros(1,2))),2));
        end
        re=(mean(loglik1)-mean(loglik2))/1000;
        
        
        
end