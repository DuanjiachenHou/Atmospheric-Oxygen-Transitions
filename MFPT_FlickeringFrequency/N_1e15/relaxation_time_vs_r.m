clear; clc; close all; 


N = 1e15;                
s = 2.03e-5;             
dr = 1.36e19;            
dd = 2.73e17;        
a1 = 0.003;           
a2 = -0.1655;           
a3 = 3.2305;            
a4 = -25.8343;         
a5 = 71.6;                


sigma_M = 1e-6;           
sigma_O = 1e-6;            
sigma = [sigma_M; sigma_O];
dt = 10;                   
T_max = 1e9;              
num_trials = 100;          
tol_ratio = 0.1;          


data = [...
    1*10^(10.65), 2.20041E+15, 1.76933E+14, 2.20041E+15, 5.60757E+18;  
    1*10^(10.70), 2.46890E+15, 1.12249E+14, 2.46890E+15, 5.20691E+18;  
    1*10^(10.80), 3.10816E+15, 6.11696E+13, 3.10816E+15, 4.45592E+18;  
    1*10^(10.90), 3.91295E+15, 3.88163E+13, 3.91295E+15, 3.76972E+18;  
    1*10^(11.00), 4.92451E+15, 2.65580E+13, 4.92611E+15, 3.14555E+18;  
    1*10^(11.10), 6.20160E+15, 1.90155E+13, 6.20160E+15, 2.58054E+18;  
    1*10^(11.20), 7.80734E+15, 1.40705E+13, 7.80736E+15, 2.07158E+18; 
    1*10^(11.30), 9.82826E+15, 1.06665E+13, 9.82888E+15, 1.61498E+18;  
    1*10^(11.40), 1.23737E+16, 8.23841E+12, 1.23738E+16, 1.20563E+18;  
    1*10^(11.50), 1.55777E+16, 6.46087E+12, 1.55777E+16, 8.34188E+17;  
    1*10^(11.60), 1.96112E+16, 5.13087E+12, 1.96112E+16, 4.66167E+17]; 

r_values = data(:,1);   
M_ss_low = data(:,2);    
O_ss_low = data(:,3);    
M_ss_high = data(:,4);   
O_ss_high = data(:,5);   

MFPT_array = zeros(length(r_values), 1); 


for i = 1:length(r_values)                          
    r_val = r_values(i);                            
    if isnan(M_ss_low(i)) || isnan(M_ss_high(i))    
        continue;                                
    end                                       
    x0 = [M_ss_low(i); O_ss_low(i)];                 
    target = [M_ss_high(i); O_ss_high(i)];          
    
    diff_state = norm(target - x0);               
    r0 = tol_ratio * diff_state;                  
    
    MFPT_trials = zeros(num_trials, 1);           
    
 
    for trial = 1:num_trials                      
        x = x0;                                 
        t = 0;                               
      
        while t < T_max                         
           
            f1 = force_MO(t, x, N, r_val, s, dr, dd, a1, a2, a3, a4, a5);  
            W = sqrt(dt) * randn(2, 1);            
            x_predict = x + f1 * dt + sigma .* W;   
         
            f2 = force_MO(t + dt, x_predict, N, r_val, s, dr, dd, a1, a2, a3, a4, a5);  
        
            x_new = x + 0.5*(f1 + f2)*dt + sigma .* W;  
            t = t + dt;                         
            x = x_new;                        
            if norm(x - target) < r0            
                break;                        
            end                            
        end                                    
        MFPT_trials(trial) = t;                
    end                                        
    
    MFPT_array(i) = mean(MFPT_trials);          
    fprintf('r = %e, MFPT = %f \n', r_val, MFPT_array(i));  
end                                            


figure;                                              
plot(log10(r_values), MFPT_array, 'o-', ...          
     'LineWidth',2, ...                             
     'DisplayName','MFPT'); hold on;                 
crit1 = 1*10^10.6127;                                
crit2 = 1*10^11.6303;                                
xline(log10(crit1), 'g--', 'LineWidth',1.5, ...        
      'DisplayName','高氧临界');                      
xline(log10(crit2), 'r--', 'LineWidth',1.5, ...      
      'DisplayName','低氧临界');                     
xlabel('log_{10}(r)', 'FontSize',12);               
ylabel('MFPT', 'FontSize',12);                      
title('不同 r 值下的平均首次穿越时间 (MFPT)', 'FontSize',14); 
legend('show');                                     
grid on;                                          


function f = force_MO(t, x, N, r, s, dr, dd, a1, a2, a3, a4, a5)  
    delta = x(2) / (dd + x(2));           
    gamma = x(2) / (dr + x(2));                   
    omega = (1 - gamma) * (1 - delta);            
    p = log10(x(2));                         
    fai = 10^(a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5);  
    f1 = 0.5 * omega * (N + r) - s * x(1) - 0.5 * fai * x(1)^0.7;  
    f2 = omega * N - (1 - omega) * r - s * x(1) - fai * x(1)^0.7;   
    f = [f1; f2];                              
end                                            
