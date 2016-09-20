
function [X1,i]=lqEstimation(X,Y,p,mu,iter,max_mu,alpha,stop)

Z=zeros(2,size(X,2));
C=zeros(2,size(X,2));
X1=X;

for i=1:iter
    dual=0;
    Z=X-Y+C/mu;
    Z=shrink(Z,mu,p);
    U=Y+Z-C/mu;
    %     [X]=rigidMotion(X,U);
    [X,Aff]=estimateA(X,U);
    for j=1:size(X,2)
        err0(j)=norm(X1(:,j)-X(:,j));
    end
    dual=max(err0);
    X1=X;
    
    P=X-Y-Z;
    C=C+mu.*P;
    if mu<max_mu 
        mu=mu.*alpha;
    end
    for j=1:size(X,2)
        err1(j)=norm(P(:,j));
    end
    primal=max(err1);
    if primal<stop && dual<stop
        break;
    end
end

end