 
 
r = [10.6, 10.65, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7];

 
lowStable = [0, 20.73, 25.58, 35.34, 43.81, 50.11, 55.74, 61.75, 66.01, 71.45, 79.18, 88.45, 100];   

 
highStable = [100, 79.27, 74.42, 64.66, 56.18, 49.89, 44.26, 38.25, 33.99, 28.54, 20.81, 11.54, 0];  

 
figure;                           

 
plot(r, lowStable, '-o', ...     
    'Color', [0.8 0.7 1], ...   
    'LineWidth', 1.5, ...        
    'MarkerSize', 6, ...         
    'DisplayName','低稳点');     
hold on;                         
 
plot(r, highStable, '-o', ...    
    'Color', [0.5 0 0.5], ...   
    'LineWidth', 1.5, ...       
    'MarkerSize', 6, ...         
    'DisplayName','高稳点');     

 
legend('Location','best');       
 
xlabel('1*10^ (r)');                     
ylabel('低稳点、高稳点 百分比(%)');  

 
title('低稳点 : 高稳点 随不同 r 的变化');  

 
grid on;                        
