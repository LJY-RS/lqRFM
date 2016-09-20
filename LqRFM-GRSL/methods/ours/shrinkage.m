function s=shrinkage(mu,n,p,s)

for i=1:3
    s=1-(p./mu).*(n.^(p-2)).*(s.^(p-1));
end
end