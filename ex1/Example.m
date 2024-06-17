syms x y nu K alpha

us1 = (exp(x)-exp(1))*cos(pi*y);
us2 = -exp(x)*sin(pi*y)/pi;
ps = 2*exp(x)*cos(pi*y);

ud1 = -(exp(x)-exp(1))*cos(pi*y);
ud2 = pi*(exp(x)-x*exp(1))*sin(pi*y);
pd = (exp(x)-x*exp(1))*cos(pi*y);

D11 = diff(us1,x,1); D12 = 1/2*(diff(us1,y,1) + diff(us2,x,1)); D22 = diff(us2,y,1);
fs1 = simplify(-2*nu*(diff(D11,x,1)+diff(D12,y,1))+diff(ps,x,1))
fs2 = simplify(-2*nu*(diff(D12,x,1)+diff(D22,y,1))+diff(ps,y,1))
gs = simplify(diff(us1,x,1)+diff(us2,y,1))

fd1 = nu*K^(-1)*ud1 + diff(pd,x,1)
fd2 = nu*K^(-1)*ud2 + diff(pd,y,1)
gd = simplify(diff(ud1,x,1)+diff(ud2,y,1))

int_ps = int(int(ps,x,0,1),y,0,1)
int_pd = int(int(pd,x,1,2),y,0,1)

eta1 = simplify(ps - 2*nu*diff(us1,x,1) - pd)
eta2 = - simplify(K^(1/2)/alpha*(diff(us1,y,1) + diff(us2,x,1)) + us2)