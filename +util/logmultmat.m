function re=logmultmat(A,B)
    nAy=size(A,1);
    nAx=size(A,2);
    nBx=size(B,2);
    for y=1:nAy
        for x=1:nBx
            re(y,x)=util.logsumexp(A(y,:)+B(:,x)',2);
        end
    end
end