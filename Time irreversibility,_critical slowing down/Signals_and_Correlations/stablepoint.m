function [x,y,t]=stablepoint(x0,y0)

% clear%清除工作空间中的所有变量
h=1;%设置时间步长h 
t=0:h:500000;%生成时间向量 t，从 0 开始，每次增加步长 h，直到 2000
N=length(t);%计算时间向量 t 的长度，赋值给 N，代表总的时间步数。
x=zeros(1,N);
y=zeros(1,N);
L=zeros(1,N);%初始化向量 x、y 和 L，大小为 1×N，元素全为 0，用于存储后续的结果。

% x(1,1) = 1.487e16;
% y(1,1) = 9.0976e+18; %给 x 和 y 的初始值赋值
x(1,1) = x0;
y(1,1) = y0; %给 x 和 y 的初始值赋值

% D = 1*10^23.4; %设置扩散系数 D
D = 1*10^20.4; %设置扩散系数 D
% a1= 0.15;
% a2= 10;
% c = 50;
% n = 0.25;
% G = 0.315;
% B = 1.4;
% a = 0.028;
% q = 0.02;
% g = 20;
% v = 0.25;
% T = 0.0001;
% k = 0.21;
% f = 0.42; %频率
% I1= 0.0275;
% I2= 0.0275;
% A = 0.025;
% b = 0.25;
% J = 2;
 


para = 10^(15.2);%设置参数 para

for i=1:N-1  %开始一个循环，从 1 到 N-1 迭代。

    [ Fx, Fy ] = force_XY(x(1,i),y(1,i),para);
    h = (10^13.3)*sqrt(Fx^2+Fy^2)^-1;
        
    x(1,i+1) = x(1,i)+(Fx)*h + sqrt(2*D*h)*randn();
    y(1,i+1) = y(1,i)+(Fy)*h + 10*sqrt(2*D*h)*randn(); %根据方程更新y，计算当前时刻y(1,i+1)的值。force_XY是一个自定义函数，返回作用于x和y的力，randn()生成正态分布的随机数模拟噪声。

    t(1,i+1) = t(1,i) + h ;
    if x(1,i+1)<0 %|| x(1,i+1)>20
        x(1,i+1)=-x(1,i+1);   %如果更新后的 x(1,i+1) 小于 0，就把它取负，保持为正值。
    end
    if y(1,i+1)<0 %|| y(1,i+1)>20
        y(1,i+1)=-y(1,i+1);
    end
end

for i=1:N  %开始一个循环，用于计算 L
    L(1,i) = sqrt(((x(1,i)-y(1,i)-20)^2)/2); %向对角线0=-x-y+20投影 %计算 L(1,i) 的值，公式为 x 和 y 相对一个常数（20）的差异平方和的平方根，表示某种投影距离。
end

lgx=log(x(1,:)+1e16);
lgy=log(y(1,:)+1e16);

figure(1)
plot(lgx,lgy);
title('Logarithmic Plot of x and y with Offset');
% set(gca,'XLim',[0,20]);
% set(gca,'YLim',[0,20]);
% 
figure(2)
plot(x,y);
title('Trajectory Plot of x versus y');
% plot(t(1,:),L(1,:));
% set(gca,'XLim',[0,250000]);
% set(gca,'YLim',[0,28.2843]);

figure(3)
plot3(x,y,-t);
title('Trajectory Plot of x versus y');

end