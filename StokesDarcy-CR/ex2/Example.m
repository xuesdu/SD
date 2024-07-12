syms x y nu K alpha_BJS

us1 = exp(x*y);
us2 = exp(x)*exp(2*y);
ps = 12*x^2*exp(y);

ud1 = exp(x*y);
ud2 = 3*x+y^2;
pd = 16*x*y^3-exp(1)-2;

D11 = diff(us1,x,1); D12 = 1/2*(diff(us1,y,1) + diff(us2,x,1)); D22 = diff(us2,y,1);
fs1 = simplify(-2*nu*(diff(D11,x,1)+diff(D12,y,1))+diff(ps,x,1))
fs2 = simplify(-2*nu*(diff(D12,x,1)+diff(D22,y,1))+diff(ps,y,1))
gs = simplify(diff(us1,x,1)+diff(us2,y,1))

% Neumann boundary condition (top and bottom)
g_b1 = nu*(diff(us1,y,1) + diff(us2,x,1))
g_b2 = 2*nu*diff(us2,y,1) - ps

fd1 = simplify(K^(-1)*ud1 + diff(pd,x,1))
fd2 = simplify(K^(-1)*ud2 + diff(pd,y,1))
gd = simplify(diff(ud1,x,1)+diff(ud2,y,1))

int_ps = int(int(ps,x,0,1),y,0,1)
int_pd = int(int(pd,x,1,2),y,0,1)

eta1 = simplify(ps - 2*nu*diff(us1,x,1) - pd)
eta2 = simplify(nu*(diff(us1,y,1) + diff(us2,x,1)) + alpha_BJS*nu^(1/2)*K^(-1/2)*us2)

us1_x = diff(us1,x,1)
us1_y = diff(us1,y,1)
us2_x = diff(us2,x,1)
us2_y = diff(us2,y,1)