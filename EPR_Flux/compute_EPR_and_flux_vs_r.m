 clear; clc; close all;                                                         

 
N     = 1*10^15;                                                               
s     = 2.03e-5;                                                               
dr    = 1.36e19;                                                             
dd    = 2.73e17;                                                             
a1    = 0.003;                                                                
a2    = -0.1655;                                                              
a3    = 3.2305;                                                                
a4    = -25.8343;                                                            
a5    = 71.6;                                                               

sigma_M = 1*10^29.4;                                                              
sigma_O = 1*10^29.4;                                                              
rel_delta = 1e-4;                                                         
eps_val   = 1.0;                                                           

 
 
 
weight_logr = [10.6, 10.65, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, ...          
               11.3, 11.4, 11.5, 11.6, 11.7];                               
p_low  = [0, 20.73, 25.58, 35.34, 43.81, 50.11, 55.74, 61.75, ...            
          66.01, 71.45, 79.18, 88.45, 100] / 100;                         
p_high = [100, 79.27, 74.42, 64.66, 56.18, 49.89, 44.26, 38.25, ...         
          33.99, 28.54, 20.81, 11.54, 0] / 100;                                

p_low_map  = containers.Map(weight_logr, p_low);                          
p_high_map = containers.Map(weight_logr, p_high);                             

 
data = [                                                                     
    1*10^9.6,    NaN,          NaN,          3.60592E+14, 1.44193E+19;        
    1*10^9.7,    NaN,          NaN,          3.90566E+14, 1.39204E+19;          
    1*10^9.8,    NaN,          NaN,          4.30739E+14, 1.33217E+19;
    1*10^9.9,    NaN,          NaN,          4.85307E+14, 1.26148E+19;
    1*10^10.0,   NaN,          NaN,          5.60103E+14, 1.17988E+19;
    1*10^10.1,   NaN,          NaN,          6.62756E+14, 1.08836E+19;
    1*10^10.2,   NaN,          NaN,          8.02910E+14, 9.89514E+18;
    1*10^10.3,   NaN,          NaN,          9.91630E+14, 8.87484E+18;
    1*10^10.4,   NaN,          NaN,          1.23966E+15, 7.86910E+18;
    1*10^10.5,   NaN,          NaN,          1.55810E+15, 6.91370E+18;
    1*10^10.6,   NaN,          NaN,          1.96114E+15, 6.02545E+18;
    1*10^10.65,  2.20041E+15,  1.76933E+14,  2.20041E+15, 5.60757E+18;
    1*10^10.7,   2.46890E+15,  1.12249E+14,  2.46890E+15, 5.20691E+18;
    1*10^10.8,   3.10816E+15,  6.11696E+13,  3.10816E+15, 4.45592E+18;
    1*10^10.9,   3.91295E+15,  3.88163E+13,  3.91295E+15, 3.76972E+18;
    1*10^11.0,   4.92451E+15,  2.65580E+13,  4.92611E+15, 3.14555E+18;
    1*10^11.1,   6.20160E+15,  1.90155E+13,  6.20160E+15, 2.58054E+18;
    1*10^11.2,   7.80734E+15,  1.40705E+13,  7.80736E+15, 2.07158E+18;
    1*10^11.3,   9.82826E+15,  1.06665E+13,  9.82888E+15, 1.61498E+18;
    1*10^11.4,   1.23737E+16,  8.23841E+12,  1.23738E+16, 1.20563E+18;
    1*10^11.5,   1.55777E+16,  6.46087E+12,  1.55777E+16, 8.34188E+17;
    1*10^11.6,   1.96112E+16,  5.13087E+12,  1.96112E+16, 4.66167E+17;
    1*10^11.7,   2.46890E+16,  4.11792E+12,  NaN,         NaN;
    1*10^11.8,   3.10816E+16,  3.33490E+12,  NaN,         NaN;
    1*10^11.9,   3.91295E+16,  2.72193E+12,  NaN,         NaN;
    1*10^12.0,   4.92611E+16,  2.23681E+12,  NaN,         NaN;
    1*10^12.1,   6.20160E+16,  1.84919E+12,  NaN,         NaN;
    1*10^12.2,   7.80736E+16,  1.53688E+12,  NaN,         NaN;
    1*10^12.3,   9.82888E+16,  1.28335E+12,  NaN,         NaN;
    1*10^12.4,   1.23738E+17,  1.07616E+12,  NaN,         NaN;
    1*10^12.5,   1.55777E+17,  9.05821E+11,  NaN,         NaN;
    1*10^12.6,   1.96112E+17,  7.65007E+11,  NaN,         NaN;
    1*10^12.7,   2.46890E+17,  6.48020E+11,  NaN,         NaN;
];                                                                       

r_values  = data(:,1);                                                       
M_ss_low  = data(:,2);                                                           
O_ss_low  = data(:,3);                                                        
M_ss_high = data(:,4);                                                         
O_ss_high = data(:,5);                                                        

 
 
 
omega_func = @(O_val) (1 - (O_val./(dr + O_val))) .* ...                  
                        (1 - (O_val./(dd + O_val)));                       
phi_func   = @(O_val) 10.^( a1*(log10(O_val)).^4 + ...                           
                           a2*(log10(O_val)).^3 + ...
                           a3*(log10(O_val)).^2 + ...
                           a4*log10(O_val) + a5 );                               

 
 
 
results = struct('r', {}, 'state', {}, 'EPR', {}, ...                       
                 'eqPoint', {}, 'V_rms', {});                                   

for i = 1:length(r_values)                                                      
    r_now = r_values(i);                                                     
    
 
    if ~isnan(M_ss_low(i)) && ~isnan(O_ss_low(i))                               
        state = 'low';                                                       
        [EPR_val, V_rms_val] = computeEPRVrms_single( ...                        
            r_now, M_ss_low(i), O_ss_low(i), ...
            N, s, sigma_M, sigma_O, omega_func, phi_func, rel_delta, eps_val);  
        
        results(end+1).r       = r_now;                                          
        results(end).state     = state;                                          
        results(end).EPR       = EPR_val;                                       
        results(end).eqPoint   = {[M_ss_low(i); O_ss_low(i)]};                 
        results(end).V_rms     = V_rms_val;                                     
    end                                                                          
    
 
    if ~isnan(M_ss_high(i)) && ~isnan(O_ss_high(i))                            
        state = 'high';                                                        
        [EPR_val, V_rms_val] = computeEPRVrms_single( ...                       
            r_now, M_ss_high(i), O_ss_high(i), ...
            N, s, sigma_M, sigma_O, omega_func, phi_func, rel_delta, eps_val);   
        
        results(end+1).r       = r_now;                                          
        results(end).state     = state;                                        
        results(end).EPR       = EPR_val;                                        
        results(end).eqPoint   = {[M_ss_high(i); O_ss_high(i)]};                
        results(end).V_rms     = V_rms_val;                                     
    end                                                                       
end                                                                             

 
 
 
results_table = struct2table(results);                                         
results_table = sortrows(results_table, 'r');                                  
results_table.logr = log10(results_table.r);                                  

 
 
 
unique_r = unique(results_table.r);                                            
num_unique = length(unique_r);                                              

EPR_weighted   = zeros(num_unique, 1);                                       
V_rms_weighted = zeros(num_unique, 1);                                       

for i = 1:num_unique                                                           
    r_current = unique_r(i);                                                
    idx = find(results_table.r == r_current);                                  
    
    if length(idx) == 2                                                       
        logr_current = log10(r_current);                                         
        
        if p_low_map.isKey(logr_current)                                         
            w_low  = p_low_map(logr_current);                                   
            w_high = p_high_map(logr_current);                                
        else                                                                   
            w_low  = 0.5;                                                       
            w_high = 0.5;                                                       
        end
        
 
        EPR_low  = results_table.EPR( idx(strcmp(results_table.state(idx),'low')) );  
        EPR_high = results_table.EPR( idx(strcmp(results_table.state(idx),'high')) );  
        
 
        V_low  = results_table.V_rms( idx(strcmp(results_table.state(idx),'low')) );  
        V_high = results_table.V_rms( idx(strcmp(results_table.state(idx),'high')) );  
        
 
        EPR_weighted(i)   = w_low * EPR_low + w_high * EPR_high;                 
        
 
        V_rms_weighted(i) = w_low * V_low  + w_high * V_high;                    
        
    else                                                                        
        EPR_weighted(i)   = results_table.EPR(idx);                             
        V_rms_weighted(i) = results_table.V_rms(idx);                         
    end                                                                       
end                                                                             


 
figure;                                                                    
crit1 = 1*10^10.6127;                                                           
crit2 = 1*10^11.6303;                                                           
plot(log10(unique_r), EPR_weighted, 'o-', 'LineWidth',2, 'DisplayName','EPR');  
hold on;                                                                     
xline(log10(crit1), 'g--', 'LineWidth',1.5, 'DisplayName','高氧临界');        
xline(log10(crit2), 'r--', 'LineWidth',1.5, 'DisplayName','低氧临界');       
xlabel('log_{10}(r)', 'FontSize',12);                                       
ylabel('EPR (加权后)', 'FontSize',12);                                         
title('不同 r 下的熵产生率（EPR）', 'FontSize',14);                     
legend('show');                                                               
grid on;                                                                      

 
 
 
V_norm = V_rms_weighted;                                                     

figure;                                                                       
crit1 = 1*10^10.6127;                                                        
crit2 = 1*10^11.6303;                                                         
plot(log10(unique_r), V_norm, 'd-', 'LineWidth',2, 'DisplayName','|V|_{rms}');  
hold on;                                                                       
xline(log10(crit1), 'g--', 'LineWidth',1.5, 'DisplayName','高氧临界');         
xline(log10(crit2), 'r--', 'LineWidth',1.5, 'DisplayName','低氧临界');         
xlabel('log_{10}(r)', 'FontSize',12);                                           
ylabel('|V|_{rms}', 'FontSize',12);                                            
title('不同 r 下的全局平均流强度 |V|_{rms}', 'FontSize',14);                  
legend('show');                                                                 
grid on;                                                                    

 

 
 
 
function [EPR_val, V_rms_val] = computeEPRVrms_single( ...
            r_val, M_ss, O_ss, N_fixed, s, sigma_M, sigma_O, ...
            omega_func, phi_func, rel_delta, eps_val)
 
    
    omega_ss = omega_func(O_ss);                                              
    phi_ss   = phi_func(O_ss);                                                 
    
 
    delta_O    = rel_delta * O_ss;                                           
    omega_plus = omega_func(O_ss + delta_O);                                 
    omega_minus= omega_func(O_ss - delta_O);                                
    phi_plus   = phi_func(O_ss + delta_O);                                  
    phi_minus  = phi_func(O_ss - delta_O);                              
    
    d_omega_dO = (omega_plus - omega_minus) / (2*delta_O);                   
    d_phi_dO   = (phi_plus - phi_minus) / (2*delta_O);                       
    
 
    dF_dM = -s - 0.35 * phi_ss * M_ss^(-0.3);                                 
    dF_dO = (0.5*(omega_plus - omega_minus)*(r_val + N_fixed) ...             
             - 0.5*(phi_plus - phi_minus)*M_ss^0.7) / (2*delta_O);           
    dG_dM = -s - 0.7 * phi_ss * M_ss^(-0.3);                                   
    dG_dO = d_omega_dO * r_val - d_omega_dO * N_fixed - d_phi_dO * M_ss^0.7;   
    
    A = [dF_dM, dF_dO; dG_dM, dG_dO];                                         
    
 
    D = diag([sigma_M^2, sigma_O^2]);                                         
    
 
    try
        Sigma = lyap(A, 2*D);                                                  
    catch
        Sigma = NaN(2);                                                        
    end
    
    if any(isnan(Sigma(:)))                                                   
        EPR_val   = NaN;                                                    
        V_rms_val = NaN;                                                     
        return;                                                             
    end
    
 
    J = A * Sigma - D;                                                       
    
 
    tmp = D \ ( J / Sigma );                                               
    EPR_matrix = J * tmp;                                                   
    EPR_val = trace(EPR_matrix) / eps_val;                                   
    
 
 
    B = J / Sigma;                                                           
    
 
    mean_V2  = trace(B' * B * Sigma);                                        
    
 
    V_rms_val = sqrt(mean_V2);                                                
end                                                                           
