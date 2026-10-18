
r = 1*10^11.6;
%N = 1*10^14.8 ;
xmin = [1e15,1e11];
xmax = [3e16,4e18];

x = linspace(xmin(1),xmax(1),501);
y = linspace(xmin(2),xmax(2),501);
[X,Y]=meshgrid(x,y);
%para = N;
para = r;
[Fx,Fy]= force_XY(X,Y,para);
v = [0,0];
contour(X,Y,Fx,v,'b');
hold on;
contour(X,Y,Fy,v,'k');
ylim([0,xmax(2)]);

dim =2;
t0 = 0; tf = 1e9;
tspan = [t0,tf];
x01 = [1e16,1e18];
x02 = [1.5e16,1.6e14];
func = @(t,x)force_MO(t,x,para);
[t1,Xt1]=ode23s(func,tspan,x01);
[t2,Xt2]=ode23s(func,tspan,x02);
X1 = Xt1(end,:);
X2 = Xt2(end,:);
plot(X1(1),X1(2),'bo');
plot(X2(1),X2(2),'bo');