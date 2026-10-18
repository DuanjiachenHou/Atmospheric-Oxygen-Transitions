clear;
clc;

logN = [14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5];     
MFPT = [20000100, 23428000, 23719300, 45081800, ...        
        61249200, 166688000, 370294000];            

nu = 1 ./ MFPT;                                            


figure;                                              
hold on;                                                     

plot(logN, nu, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);     

xlabel('log_{10}(N)', 'FontSize', 12);                      
ylabel('Flickering frequency \nu', 'FontSize', 12);         

set(gca, 'FontSize', 12, 'LineWidth', 1.2);                  
xlim([14.89, 15.61]);                                     
box on;                                                     
% set(gca, 'YScale', 'log');                              


x_left  = 14.9;                                          
x_right = 15.6;                                            

xline(x_left,  'r--', 'LineWidth', 1.5);                      
xline(x_right, 'g--', 'LineWidth', 1.5);                  

title('Flickering frequency vs log_{10}(N)', 'FontSize', 13); 

hold off;                                                    
