function[r] = hodnost(matica, a)
[s, v, d] = svd(matica);
vl_hodnoty = diag(v);
r = max(find(vl_hodnoty>a));
end