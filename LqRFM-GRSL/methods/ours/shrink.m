function Z=shrink(Z,u,p)

Aa=(2./u.*(1-p)).^(1/(2-p));
ha=Aa+p./u.*(Aa.^(p-1));

for i=1:size(Z,2)
    n=norm(Z(:,i));
    w=0;
    if n>ha
        w = shrinkage(u, n, p, (Aa/n + 1.0)/2.0);
    end
    Z(:,i)=Z(:,i)*w;
end

end