clear; clc; close all;  


sigma_M = 1e-6;         
sigma_O = 1e-6;           
sigma = [sigma_M; sigma_O];


dt = 10;              
T_max = 1e8;              
num_trials = 100;        

tol_ratio = 0.1;          

data = [...
    1*10^15.0, 1.47783E+16, 6.82215E+12,    1.47783E+16,    9.16684E+17;   
    1*10^15.1, 1.47781E+16, 9.74415E+12,    1.47783E+16,    1.46898E+18;  
    1*10^15.2, 1.47781E+16, 1.44235E+13,    1.47783E+16,    2.11124E+18;   
    1*10^15.3, 1.47782E+16, 2.24673E+13,    1.47783E+16,    2.86041E+18;   
    1*10^15.4, 1.47784E+16, 3.79197E+13,    1.47783E+16,    3.72684E+18;   
    1*10^15.5, 1.47784E+16, 7.46981E+13,    1.47783E+16,    4.71910E+18;   
    1*10^15.6, 1.47785E+16, 2.86412E+14,    1.47783E+16,    5.84536E+18];  

 
N_values  = data(:,1);    
M_ss_low  = data(:,2);    
O_ss_low  = data(:,3);    
M_ss_high = data(:,4);    
O_ss_high = data(:,5);    

 
MFPT_array = zeros(length(N_values),1);   

 
for i = 1:length(N_values)        
    N = N_values(i);             
    x0 = [M_ss_low(i); O_ss_low(i)];     
    target = [M_ss_high(i); O_ss_high(i)]; 
    
    
    diff_state = norm(target - x0);  
    r0 = tol_ratio * diff_state;     
    
    MFPT_trials = zeros(num_trials,1);  
    
     
    for trial = 1:num_trials
        x = x0;           
        t = 0;             
        
         
        while t < T_max
            
            f1 = force_MO(t, x, N);     
             
            W = sqrt(dt) * randn(2,1);     
           
            x_predict = x + f1*dt + sigma.*W;  
           
            f2 = force_MO(t + dt, x_predict, N); 
           
            x_new = x + 0.5*(f1 + f2)*dt + sigma.*W;  
            t = t + dt;                 
            x = x_new;                
            
           
            if norm(x - target) < r0
                break;                
            end
        end
        MFPT_trials(trial) = t;  
    end
    
 
    MFPT_array(i) = mean(MFPT_trials);  
    fprintf('N = %e, MFPT = %f \n', N, MFPT_array(i));  
end



figure;                                  
crit1 = 1*10^14.8999;                     
crit2 = 1*10^15.6026;                     
plot(log10(N_values), MFPT_array, 'o-', ...  
    'LineWidth',2, 'DisplayName','MFPT');  
hold on;                                   
xline(log10(crit1), 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');              
xline(log10(crit2), 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');             
xlabel('log_{10}(N)', 'FontSize',12);     
ylabel('MFPT', 'FontSize',12);            
title('不同 N 下的平均首次穿越时间 (MFPT)', 'FontSize',14);  
legend('show');                           
grid on;                                 


function f = force_MO(t, x, N)
    
    r = 3e11;          
    s = 2.03e-5;       
    dr = 1.36e19;      
    dd = 2.73e17;      
    a1 = 0.003;       
    a2 = -0.1655;      
    a3 = 3.2305;       
    a4 = -25.8343;     
    a5 = 71.6;         
    
  
    delta = x(2) / (dd + x(2));  
    gamma = x(2) / (dr + x(2));    
    omega = (1 - gamma) * (1 - delta);  
    p = log10(x(2));            
    
    fai = 10^(a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5);  
   
    f1 = 0.5 * omega * (N + r) - s * x(1) - 0.5 * fai * x(1)^0.7; 
    f2 = omega * N - (1 - omega) * r - s * x(1) - fai * x(1)^0.7;  
    f = [f1; f2];  
end
