function [Aff,idxF]=lqRFM(X,Y,scores,t1,t2)

p=0.2;
mu=0.0003;
alpha=1.65;
Nmax=100;
Nmin=50;
iter=30;
max_mu=1e6;
stop=1e-8;

[Xn, Yn, ~]=norm2(X,Y);
nsample=max(min(size(Xn,1)*0.3,Nmax),Nmin);
[vals,indx] = sort(scores);
xsample=Xn(indx(1:nsample),:)';
ysample=Yn(indx(1:nsample),:)';

Xsample=X(indx(1:nsample),:)';
Ysample=Y(indx(1:nsample),:)';

[x1_,ite]=lqEstimation(ysample,xsample,p,mu,iter,max_mu,alpha,stop);
err=sum((x1_-xsample).^2);
inliers=err<t1*t1;

[~,idx]=find(inliers==1);
X1=Xsample(:,idx);
X2=Ysample(:,idx);
[Xe,Aff]=estimateA(X2,X1);
vv=(Xe-X1);
m0 = sqrt(sum(vv(1,:).^2+vv(2,:).^2)/(size(vv,2)));
if m0>0.5&m0<5
    tt=3*m0;
else
    tt=t2;
end

X_=[Aff(1:2)';Aff(4:5)']*Y'+repmat([Aff(3);Aff(6)],1,size(Y',2));
E=X_-X';
dstE = sqrt(E(1,:).^2+E(2,:).^2);
inliersF=dstE<tt;
[~,idxF]=find(inliersF==1);

X1=X(idxF,:)';
X2=Y(idxF,:)';
[~,Aff]=estimateA(X2,X1);

X_=[Aff(1:2)';Aff(4:5)']*Y'+repmat([Aff(3);Aff(6)],1,size(Y',2));
E=X_-X';
dstE = sqrt(E(1,:).^2+E(2,:).^2);
inliersF=dstE<t2;
[~,idxF]=find(inliersF==1);