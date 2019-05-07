function [d] = kl_normaltrunc (mu1,sig1,mu2,sig2,a,b)
    %norm1=sum(normpdf([a:b],mu1,sig1));
    %norm2=sum(normpdf([a:b],mu2,sig2));
    syms x
    p=1/(sig1*sqrt(2*pi))*exp(-1/2*((x-mu1)/sig1)^2);
    q=1/(sig2*sqrt(2*pi))*exp(-1/2*((x-mu2)/sig2)^2);
    f=p*log(p/q);
    re=int(f,a,b);
    d=eval(re);
end