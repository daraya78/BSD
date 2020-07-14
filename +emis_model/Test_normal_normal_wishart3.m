
clear all
mu=ones(1,2);
cova=[10 4;4 5];
ndata=1000;

data=mvnrnd(mu,cova,ndata);
plot(data(:,1),data(:,2),'x')

a0=1;
b0=1;
m0=zeros(1,2);
v0=2;
W0=eye(2);

N=ndata;
mk=rand(1,2);
ak=a0;
bk=b0;
gamma=ones(ndata,1);
xbar=mean(data);
XC = bsxfun(@minus,data,xbar);
Sk = bsxfun(@times, XC, gamma)'*XC / N;
for j=1:ndata
    vk=N+v0;
    Wk=inv((W0^-1)+N*Sk);
    Rk=diag(ak.*bk)+ndata*(vk*Wk);
    mk=(m0*diag(ak.*bk)+xbar*N*(vk*Wk))*(Rk)^-1;
    ak=a0+1/2;
    bk=(1./b0+1/2*((mk-m0)*(mk-m0)'+trace(inv(Rk)))).^-1;
    bkk(j)=bk;
end

plot(data(:,1),data(:,2),'.')
hold on
h1=MyEllipse(Rk^-1*4,mk);
h1.LineWidth=2;
hold on
h2=MyEllipse((Wk*vk)^-1*4,mk);
h2.LineWidth=2;
h3=MyEllipse((ak*bk*eye(2))^-1*4,mk);
h3.LineWidth=2;
axis([-10 10 -10 10])

bk











