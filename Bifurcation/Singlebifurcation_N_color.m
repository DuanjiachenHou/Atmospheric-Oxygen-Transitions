 
clear;                                 
clc;                                    

 
r = 3e11;                               
s = 2.03e-5;                            
dr = 1.36e19;                           
dd = 2.73e17;                         
a1 = 0.003;                            
a2 = -0.1655;                           
a3 = 3.2305;                           
a4 = -25.8343;                     
a5 = 71.6;                           

 
Y_range = linspace(10,20,100);          

 
Y = Y_range;                           
delta = 10.^Y ./ (dd + 10.^Y);       
gamma = 10.^Y ./ (dr + 10.^Y);        
omega = (1 - gamma) .* (1 - delta);   
p = Y;                                 
ap = a1.*p.^4 + a2.*p.^3 + a3.*p.^2 + a4.*p + a5;   
fai = 10.^ap;                          
N = fai ./ (omega .* (s / r).^(0.7));  

 
x = log10(N);                         

 
 
idx1 = (x <= 15.6051) & (Y < 14.6465);  

 
idx2 = (x >= 14.8929) & (x <= 15.6051) & (Y >= 14.6465) & (Y <= 17.4747);

 
idx3 = (x > 14.8929) & (Y > 17.4747);

 
figure;                              
hold on;                              

 
plot(x(idx1), Y(idx1), 'm-', 'LineWidth',1.5);   

 
plot(x(idx2), Y(idx2), 'k--','LineWidth',1.5);  

 
plot(x(idx3), Y(idx3), 'g-', 'LineWidth',1.5);   

hold off;                          

 
 
title('Bifurcation structure under variation of N');      
xlabel('log_{10}(N)');                 % X 轴标签
ylabel('log_{10}(O)');                           % Y 轴标签
grid on;                               % 打开网格

 
 


legend({'Low-oxygen state', 'Unstable branch', 'High-oxygen state'}, ...
       'Location', 'Best');           %
