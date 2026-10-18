function f = force_MO(t,x,para)
N = para;
% N = 3.75e15;
% N = 1e15;
r = 3e11;
s = 2.03e-5;
% w = 6e-9;
dr = 1.36e19;
dd = 2.73e17;
a1 = 0.003; a2=-0.1655; a3=3.2305; a4=-25.8343; %a5=71.5398;
a5=71.6;
delta = x(2)/(dd+x(2));
gamma = x(2)/(dr+x(2));
omega = (1-gamma)*(1-delta);
p = log10(x(2));
fai = 10^(a1*p^4+a2*p^3+a3*p^2+a4*p+a5);
f = [0.5*omega*(N+r)-s*x(1)-0.5*fai*x(1)^0.7;
    omega*N-(1-omega)*r-s*x(1)-fai*x(1)^0.7];
