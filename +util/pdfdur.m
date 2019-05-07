 function re = pdfdur(emision,durmax)
    T=size(emision,1);
    re=zeros(T,durmax);
    re(1,1)=emision(1);
    re(2,:)=0;
    for t=2:T
        
        re(t,1)=emision(t);
        for d=2:durmax
            re(t,d)=re(t-1,d-1)+re(t,1);
        end
    end
 end