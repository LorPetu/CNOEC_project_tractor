ns=100;
nu=2;
np=ceil(ns/nu);
u_in=ones(np,1);
test=zeros(ns,2);
for ind=2:ns+1
    test(ind,2)                 =   u_in(floor(ind*ns/Nu));
    test(ind,1)                 =   ind

end











