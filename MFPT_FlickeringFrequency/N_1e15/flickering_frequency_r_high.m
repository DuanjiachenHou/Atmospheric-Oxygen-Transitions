clear;
clc;

logr = [10.65, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6];
MFPT = [56289500, 40202400, 29005500, 35711300, 28499400, ...             
        22523700, 14365400, 7579200, 5731400, 4083850, 2207900];         


nu = 1 ./ MFPT;                                                               


figure;                                                                     
hold on;                                                                     

plot(logr, nu, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);                  

xlabel('log_{10}(r)', 'FontSize', 12);                                      
ylabel('Flickering frequency \nu', 'FontSize', 12);                        
title('Flickering frequency vs r (\it含临界线)', 'FontSize', 13);           

set(gca, 'FontSize', 12, 'LineWidth', 1.2);                                  
xlim([10.6, 11.6]);                                                        
box on;                                                                     



x_left = 10.65;                                                           
x_right = 11.6;                                                             

xline(x_left, 'g--', 'LineWidth', 1.5);                                  
xline(x_right, 'r--', 'LineWidth', 1.5);                                      

hold off;                                                                   