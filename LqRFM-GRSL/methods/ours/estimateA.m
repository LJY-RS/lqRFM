function [X1,Aff]=estimateA(X1,X2)

n=size(X1,2);
nrows_M=2*n;
ncols_M=6;
M=zeros(nrows_M,ncols_M);
B=zeros(nrows_M,1);
for i=1:n
    x1=X2(1,i);
    y1=X2(2,i);
    x2=X1(1,i);
    y2=X1(2,i);

    M_=[x2, y2, 1, 0, 0, 0;
        0, 0, 0, x2, y2, 1];
    
    b=[x1;y1];
    row_ini=i*2-1;
    row_end=i*2;
        
    M(row_ini:row_end,:)=M_;
    B(row_ini:row_end,1)=b;
end

Aff=inv(M'*M)*(M'*B);
X1=[Aff(1:2)';Aff(4:5)']*X1+repmat([Aff(3);Aff(6)],1,size(X1,2));

