 

clear; clc; close all;                                                        

 
r  = 3*10^11;        
s  = 2.03e-5;        
dr = 1.36e19;       
dd = 2.73e17;       
a1 = 0.003;          
a2 = -0.1655;        
a3 = 3.2305;         
a4 = -25.8343;       
a5 = 71.6;            

 
sigma_M = 1*10^29.4;       
sigma_O = 1*10^29.4;       

eps_val = 1.0;        

 
weight_logN = [14.9, 15.0, 15.1, 15.2, 15.3, 15.4, 15.5, 15.6];               
w1 = [93.10, 77.71, 65.81, 60.37, 53.29, 44.97, 32.18, 9.28];                  
w2 = [6.90, 22.27, 34.18, 39.63, 46.67, 55.01, 67.78, 90.72];                 
p_low_weights  = w1 / 100;                                                    
p_high_weights = w2 / 100;                                                     

 
 
data = [                                                                       
    1*10^14.0, 1.47783E+16, 4.66623E+11,    NaN,            NaN;
    1*10^14.1, 1.47783E+16, 5.86664E+11,    NaN,            NaN;
    1*10^14.2, 1.47783E+16, 7.41629E+11,    NaN,            NaN;
    1*10^14.3, 1.47783E+16, 9.43449E+11,    NaN,            NaN;
    1*10^14.4, 1.47783E+16, 1.20895E+12,    NaN,            NaN;
    1*10^14.5, 1.47783E+16, 1.56228E+12,    NaN,            NaN;
    1*10^14.6, 1.47783E+16, 2.03885E+12,    NaN,            NaN;
    1*10^14.7, 1.47783E+16, 2.69193E+12,    NaN,            NaN;
    1*10^14.8, 1.47811E+16, 3.60494E+12,    NaN,            NaN;
    1*10^14.9, 1.47784E+16, 4.90558E+12,    1.47783E+16,    3.72183E+17;
    1*10^15.0, 1.47783E+16, 6.82215E+12,    1.47783E+16,    9.16684E+17;
    1*10^15.1, 1.47781E+16, 9.74415E+12,    1.47783E+16,    1.46898E+18;
    1*10^15.2, 1.47781E+16, 1.44235E+13,    1.47783E+16,    2.11124E+18;
    1*10^15.3, 1.47782E+16, 2.24673E+13,    1.47783E+16,    2.86041E+18;
    1*10^15.4, 1.47784E+16, 3.79197E+13,    1.47783E+16,    3.72684E+18;
    1*10^15.5, 1.47784E+16, 7.46981E+13,    1.47783E+16,    4.71910E+18;
    1*10^15.6, 1.47785E+16, 2.86412E+14,    1.47783E+16,    5.84536E+18;
    1*10^15.7,       NaN,         NaN,      1.47783E+16,    7.11396E+18;
    1*10^15.8,       NaN,         NaN,      1.47783E+16,    8.53375E+18;
    1*10^15.9,       NaN,         NaN,      1.47783E+16,    1.01143E+19;
    1*10^16.0,       NaN,         NaN,      1.47783E+16,    1.18657E+19;
    1*10^16.1,       NaN,         NaN,      1.47783E+16,    1.37994E+19;
    1*10^16.2,       NaN,         NaN,      1.47783E+16,    1.59272E+19;
    1*10^16.3,       NaN,         NaN,      1.47783E+16,    1.82621E+19;
    1*10^16.4,       NaN,         NaN,      1.47783E+16,    2.08180E+19;
    1*10^16.5,       NaN,         NaN,      1.47783E+16,    2.36098E+19;
];                                                                             

 
N_values  = data(:,1);    
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

 
results = struct('N', {}, 'state', {}, 'EPR', {}, 'EPR_normalized', {}, ...    
                 'Sigma', {}, 'eqPoint', {}, 'V_rms', {});                     

relative_delta = 1e-4;                                                         
num_N = length(N_values);                                                      

 
for j = 1:num_N                                                              
    if ~isnan(M_ss_low(j)) && ~isnan(O_ss_low(j))                             
        state = 'low';                                                         
        N = N_values(j);                                                      
        M_ss = M_ss_low(j);                                                   
        O_ss = O_ss_low(j);                                                   
        
        fprintf('\n--- 计算 N = 1*10^%.1f, 稳态点: %s ---\n', log10(N), state);
        fprintf('M_ss = %.6e, O_ss = %.6e\n', M_ss, O_ss);                    
        
        omega_ss = omega_func(O_ss);                                          
        phi_ss   = phi_func(O_ss);                                            
        fprintf('omega_ss = %.6e, phi_ss = %.6e\n', omega_ss, phi_ss);         
        
 
        dF_dM = -s - 0.35*phi_ss*M_ss^(-0.3);                                 
        
        delta_O = relative_delta * O_ss;                                       
        omega_plus  = omega_func(O_ss + delta_O);                              
        omega_minus = omega_func(O_ss - delta_O);                              
        phi_plus    = phi_func(O_ss + delta_O);                                
        phi_minus   = phi_func(O_ss - delta_O);                                
        
        dF_dO = (0.5*(omega_plus - omega_minus)*(N + r) - ...                  
                 0.5*(phi_plus - phi_minus)*M_ss^0.7) / (2*delta_O);           
        
        dG_dM = -s - 0.7*phi_ss*M_ss^(-0.3);                                   
        
        d_omega_dO = (omega_plus - omega_minus)/(2*delta_O);                   
        d_phi_dO   = (phi_plus  - phi_minus )/(2*delta_O);                   
        
        dG_dO = d_omega_dO*N - d_omega_dO*r - d_phi_dO*M_ss^0.7;            
        
        Ai = [dF_dM, dF_dO; dG_dM, dG_dO];                                    
        disp('雅可比矩阵 Ai:'); disp(Ai);                                    
        
 
        Di = diag([sigma_M^2, sigma_O^2]);                                     
        disp('噪声矩阵 Di:'); disp(Di);                                      
        
 
        eig_Ai = eig(Ai);                                                   
        disp('特征值:'); disp(eig_Ai);                                        
        if any(real(eig_Ai) >= 0)                                             
            fprintf('警告：雅可比矩阵不满足 Hurwitz 条件！\n');              
        end                                                                   
        
 
        try                                                                    
            Sigma = lyap(Ai, 2*Di);                                         
        catch ME                                                             
            fprintf('Lyapunov 方程求解失败: %s\n', ME.message);              
            Sigma = NaN(2);                                               
        end                                                                 
        
        if ~any(isnan(Sigma(:)))                                             
            disp('协方差矩阵 Σ:'); disp(Sigma);                               
            
 
            J = Ai * Sigma - Di;                                             
            
 
 
            tmp = Di \ ( J / Sigma );                                       
            EPR_matrix = J * tmp;                                          
            EPR = trace(EPR_matrix) / eps_val;                        
            
            fprintf('EPR = %.6e\n', EPR);                                 
            
            EPR_normalized = EPR / N;                                       
            fprintf('归一化 EPR = %.6e\n', EPR_normalized);               
            
 
 
            B = J / Sigma;                                                   
            
 
            mean_V2 = trace(B' * B * Sigma);                                
            
            V_rms = sqrt(mean_V2);                                        
            fprintf('平均流强度 |V|_rms = %.6e\n', V_rms);                 
            
 
            results(end+1).N = N;                                         
            results(end).state = state;                                     
            results(end).EPR = EPR;                                      
            results(end).EPR_normalized = EPR_normalized;                   
            results(end).Sigma = {Sigma};                                    
            results(end).eqPoint = { [M_ss; O_ss] };                         
            results(end).V_rms = V_rms;                                      
        else                                                                
            fprintf('低稳态 N = 1*10^%.1f 的 Sigma 含 NaN，跳过。\n', log10(N));
            results(end+1).N = N;                                             
            results(end).state = state;                                     
            results(end).EPR = NaN;                                         
            results(end).EPR_normalized = NaN;                              
            results(end).Sigma = {NaN(2)};                              
            results(end).eqPoint = { [M_ss; O_ss] };                        
            results(end).V_rms = NaN;                                      
        end                                                              
    end                                                                      
end                                                                          

 
for j = 1:num_N                                                                
    if ~isnan(M_ss_high(j)) && ~isnan(O_ss_high(j))                           
        state = 'high';                                                       
        N = N_values(j);                                                      
        M_ss = M_ss_high(j);                                                 
        O_ss = O_ss_high(j);                                                 
        
        fprintf('\n--- 计算 N = 1*10^%.1f, 稳态点: %s ---\n', log10(N), state);
        fprintf('M_ss = %.6e, O_ss = %.6e\n', M_ss, O_ss);                     
        
        omega_ss = omega_func(O_ss);                                          
        phi_ss   = phi_func(O_ss);                                          
        fprintf('omega_ss = %.6e, phi_ss = %.6e\n', omega_ss, phi_ss);        
        
 
        dF_dM = -s - 0.35*phi_ss*M_ss^(-0.3);                               
        delta_O = relative_delta * O_ss;                                      
        omega_plus  = omega_func(O_ss + delta_O);                        
        omega_minus = omega_func(O_ss - delta_O);                            
        phi_plus    = phi_func(O_ss + delta_O);                            
        phi_minus   = phi_func(O_ss - delta_O);                               
        
        dF_dO = (0.5*(omega_plus - omega_minus)*(N + r) - ...                 
                 0.5*(phi_plus - phi_minus)*M_ss^0.7) / (2*delta_O);     
        
        dG_dM = -s - 0.7*phi_ss*M_ss^(-0.3);                                  
        d_omega_dO = (omega_plus - omega_minus)/(2*delta_O);                  
        d_phi_dO   = (phi_plus  - phi_minus )/(2*delta_O);                   
        dG_dO = d_omega_dO*N - d_omega_dO*r - d_phi_dO*M_ss^0.7;            
        
        Ai = [dF_dM, dF_dO; dG_dM, dG_dO];                                  
        disp('雅可比矩阵 Ai:'); disp(Ai);                                   
        
 
        Di = diag([sigma_M^2, sigma_O^2]);                                  
        disp('噪声矩阵 Di:'); disp(Di);                                     
        
 
        eig_Ai = eig(Ai);                                                 
        disp('特征值:'); disp(eig_Ai);                                      
        if any(real(eig_Ai) >= 0)                                            
            fprintf('警告：雅可比矩阵不满足 Hurwitz 条件！\n');             
        end                                                                  
        
 
        try                                                                   
            Sigma = lyap(Ai, 2*Di);                                           
        catch ME                                                           
            fprintf('Lyapunov 方程求解失败: %s\n', ME.message);               
            Sigma = NaN(2);                                                
        end                                                                  
        
        if ~any(isnan(Sigma(:)))                                               
            disp('协方差矩阵 Σ:'); disp(Sigma);                                
            
 
            J = Ai * Sigma - Di;                                            
            
 
            tmp = Di \ ( J / Sigma );                                      
            EPR_matrix = J * tmp;                                           
            EPR = trace(EPR_matrix) / eps_val;                              
            
            fprintf('EPR = %.6e\n', EPR);                                    
            
            EPR_normalized = EPR / N;                                        
            fprintf('归一化 EPR = %.6e\n', EPR_normalized);                  
            
 
            B = J / Sigma;                                                 
            mean_V2 = trace(B' * B * Sigma);                                
            V_rms = sqrt(mean_V2);                                            
            fprintf('平均流强度 |V|_rms = %.6e\n', V_rms);                   
            
 
            results(end+1).N = N;                                           
            results(end).state = state;                                     
            results(end).EPR = EPR;                                      
            results(end).EPR_normalized = EPR_normalized;                   
            results(end).Sigma = {Sigma};                                    
            results(end).eqPoint = { [M_ss; O_ss] };                        
            results(end).V_rms = V_rms;                                      
        else                                                                  
            fprintf('高稳态 N = 1*10^%.1f 的 Sigma 含 NaN，跳过。\n', log10(N));
            results(end+1).N = N;                                            
            results(end).state = state;                             
            results(end).EPR = NaN;                                        
            results(end).EPR_normalized = NaN;                               
            results(end).Sigma = {NaN(2)};                                 
            results(end).eqPoint = { [M_ss; O_ss] };                         
            results(end).V_rms = NaN;                                        
        end                                                                    
    end                                                                      
end                                                                

 
results_table = struct2table(results);                                        
if any(results_table.N <= 0)                                               
    error('N 值中存在非正数，请检查数据。');                                 
end                                                                     
sorted_table = sortrows(results_table, 'N');                                 
sorted_table.logN = log10(sorted_table.N);                                  

 
unique_N = unique(sorted_table.N);                                          
num_unique = length(unique_N);                                               

EPR_weighted   = zeros(num_unique,1);                                
V_rms_weighted = zeros(num_unique,1);                                        

for i = 1:num_unique                                                         
    current_N = unique_N(i);                                                 
    idx = find(sorted_table.N == current_N);                                 
    
    if length(idx) == 2                                                     
        current_logN = log10(current_N);                                      
        
 
        if ismember(current_logN, weight_logN)                               
            weight_idx = find(weight_logN == current_logN);                   
            w_low  = p_low_weights(weight_idx);                           
            w_high = p_high_weights(weight_idx);                             
        else                                                                  
            w_low  = 0.5;                                                     
            w_high = 0.5;                                                    
        end                                                                  
        
 
        idx_low  = idx(strcmp(sorted_table.state(idx), 'low'));              
        idx_high = idx(strcmp(sorted_table.state(idx), 'high'));              
        
        EPR_low  = sorted_table.EPR(idx_low);                                 
        EPR_high = sorted_table.EPR(idx_high);                                 
        
        V_low  = sorted_table.V_rms(idx_low);                                
        V_high = sorted_table.V_rms(idx_high);                                
        
 
        EPR_weighted(i) = w_low * EPR_low + w_high * EPR_high;                
        
 
        V_rms_weighted(i) = w_low * V_low + w_high * V_high;                 
        
    else                                                                    
        EPR_weighted(i)   = sorted_table.EPR(idx);                           
        V_rms_weighted(i) = sorted_table.V_rms(idx);                          
    end                                                                        
end                                                                           

 
figure;                                                                     
crit1 = 1*10^14.8999;                                                          
crit2 = 1*10^15.6026;                                                         
plot(log10(unique_N), EPR_weighted, 'o-', 'LineWidth', 2);                    
xline(log10(crit1), 'r--', 'LineWidth',1.5, 'DisplayName','低氧临界');         
xline(log10(crit2), 'g--', 'LineWidth',1.5, 'DisplayName','高氧临界');         
xlabel('log_{10}(N)', 'FontSize', 12);                                        
ylabel('EPR', 'FontSize', 12);                                               
title('不同 N 下的熵产生率（EPR）', 'FontSize', 14);                        
grid on;                                                                     

 
V_norm = V_rms_weighted;                                                     

figure;                                                                       
crit1 = 1*10^14.8999;                                                        
crit2 = 1*10^15.6026;                                                   
plot(log10(unique_N), V_norm, 'd-', ...                                    
     'LineWidth',2, 'DisplayName','|V|_{rms}'); hold on;                     
xline(log10(crit1), 'r--', 'LineWidth',1.5, ...                             
      'DisplayName','低氧临界');                                              
xline(log10(crit2), 'g--', 'LineWidth',1.5, ...                              
      'DisplayName','高氧临界');                                              
xlabel('log_{10}(N)', 'FontSize',12);                                         
ylabel('|V|_{rms}', 'FontSize',12);                                           
title('不同 N 下的平均流强度 |V|_{rms}', 'FontSize',14);                       
legend('show');                                                                
grid on;                                                                   

