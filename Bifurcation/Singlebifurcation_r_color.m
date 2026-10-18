 
clear;                       
clc;                         
 

 
N = 1e15;                     
s = 2.03e-5;                  
dr = 1.36e19;                 
dd = 2.73e17;                 
a1 = 0.003;                    
a2 = -0.1655;                 
a3 = 3.2305;                  
a4 = -25.8343;                
a5 = 71.6;                    

 
Y_range = linspace(10,20,1000);   

 
 
Y = Y_range;                                          
delta = 10.^Y ./ (dd + 10.^Y);                          
gamma = 10.^Y ./ (dr + 10.^Y);                       
omega = (1 - gamma) .* (1 - delta);                    
p = Y;                                                
ap = a1.*p.^4 + a2.*p.^3 + a3.*p.^2 + a4.*p + a5;        
fai = 10.^ap;                                          
r = s ./ ((fai ./ (omega * N)).^(1/0.7));              
x = log10(r);                                           

 
 
idx1 = (x <= 11.6303) & (Y > 17.4747);

 
idx2 = (x >= 10.6127) & (x <= 11.6303) & (Y >= 14.6465) & (Y <= 17.4747);

 
idx3 = (x > 10.6127) & (Y < 14.6465);

 
figure;                      
hold on;                     

 
plot(x(idx1), Y(idx1), 'g-', 'LineWidth', 1.5);   

 
plot(x(idx2), Y(idx2), 'k--', 'LineWidth', 1.5); 

 
plot(x(idx3), Y(idx3), 'm-', 'LineWidth', 1.5);   

hold off;                     

 
 
title('Bifurcation structure under variation of r'); 
xlabel('log_{10}(r)');        
ylabel('log_{10}(O)');                 
grid on;                      

 

legend({'Low-oxygen state', 'Unstable branch', 'High-oxygen state'}, ...
       'Location', 'Best');           %
