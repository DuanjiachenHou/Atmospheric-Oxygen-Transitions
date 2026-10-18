 clear; clc; close all;                                                       

 
sigma_M = 1e14;                                                               
sigma_O = 1e14;                                                                 
sigma   = [sigma_M; sigma_O];                                                  

 
dt        = 1e2;                                                              
T_max     = 1e8;                                                              
num_trials = 100;                                                                
burn_in_ratio = 0.02;                                                        

 
stay_in_ratio_high  = 0.25;                                                     
stay_out_ratio_high = 0.35;                                                     
assert(stay_out_ratio_high>stay_in_ratio_high, '外圈必须大于内圈');                

 
cross_ratio = 0.55;                                                             
cone_cos_min = 0.0;                                                         

 
use_buffer_band   = true;                                                    
buffer_band_ratio = 0.12;                                                    

 
progress_print_every = 0.05;                                                   
max_traj_points = 5000;                                                      

 
crit1 = 1*10^14.8999;                                                         
crit2 = 1*10^15.6026;                                                          

rng('shuffle');                                                               

 
 
data = [ ...                                                                 
    1*10^14.9, 1.47784E+16, 4.90558E+12, 1.47783E+16, 3.72183E+17, NaN, NaN;   
    1*10^15.0, 1.47783E+16, 6.82215E+12, 1.47783E+16, 9.16684E+17, NaN, NaN;  
    1*10^15.1, 1.47781E+16, 9.74415E+12, 1.47783E+16, 1.46898E+18, NaN, NaN;   
    1*10^15.2, 1.47781E+16, 1.44235E+13, 1.47783E+16, 2.11124E+18, NaN, NaN;   
    1*10^15.3, 1.47782E+16, 2.24673E+13, 1.47783E+16, 2.86041E+18, NaN, NaN;   
    1*10^15.4, 1.47784E+16, 3.79197E+13, 1.47783E+16, 3.72684E+18, NaN, NaN;   
    1*10^15.5, 1.47784E+16, 7.46981E+13, 1.47783E+16, 4.71910E+18, NaN, NaN;   
    1*10^15.6, 1.47785E+16, 2.86412E+14, 1.47783E+16, 5.84536E+18, NaN, NaN]; 
 
 
N_values  = data(:,1);                                                       
M_low     = data(:,2);                                                       
O_low     = data(:,3);                                                       
M_high    = data(:,4);                                                        
O_high    = data(:,5);                                                        
M_sad_raw = data(:,6);                                                        
O_sad_raw = data(:,7);                                                      

 
FPT_down_all   = NaN(length(N_values),1);                                      
hit_flag_all   = false(length(N_values),1);                                    
frac_high_all  = zeros(length(N_values),1);                                   
frac_low_near_all  = zeros(length(N_values),1);                               
frac_high_near_all = zeros(length(N_values),1);                               
frac_buffer_all    = zeros(length(N_values),1);                               

 
for i = 1:length(N_values)                                                   
    N = N_values(i);                                                           
    high_pt = [M_high(i); O_high(i)];                                        
    low_pt  = [M_low(i);  O_low(i)];                                         

 
    if isnan(M_sad_raw(i)) || isnan(O_sad_raw(i))                              
        sad_pt = 0.5*(high_pt + low_pt);                                       
    else                                                                       
        sad_pt = [M_sad_raw(i); O_sad_raw(i)];                                
    end

 
    D_HS = norm(sad_pt - high_pt);                                            
    D_HL = norm(low_pt  - high_pt);                                           

 
    r_in_high  = stay_in_ratio_high  * D_HS;                                 
    r_out_high = stay_out_ratio_high * D_HS;                                  
    r_cross    = cross_ratio          * D_HS;                                 

 
    buf_w = buffer_band_ratio * D_HL;                                       

 
    FPT_trials   = NaN(num_trials,1);                                        
    hit_trials   = false(num_trials,1);                                      
    near_low_trials   = zeros(num_trials,1);                                
    near_high_trials  = zeros(num_trials,1);                               
    buffer_trials     = zeros(num_trials,1);                                 
    high_stay_trials  = zeros(num_trials,1);                                 
    total_time_trials = zeros(num_trials,1);                                 

 
    steps_total = floor(T_max/dt);                                         
    sample_step = max(1, floor(steps_total/max_traj_points));               

    for trial = 1:num_trials                                                 
        t = 0;                                                              
        x = high_pt;                                                        
        burn_steps = floor(burn_in_ratio*steps_total);                      
        next_print = progress_print_every;                                 

 
        if trial==1                                                          
            Xs = nan(1, ceil(steps_total/sample_step));                       
            Os = nan(1, ceil(steps_total/sample_step));                       
            si = 0;                                                         
        end

 
        in_high = true;                                                  
        t_high    = 0;                                                      
        t_near_L  = 0;                                                     
        t_near_H  = 0;                                                   
        t_buffer  = 0;                                                     
        t_alarm   = NaN;                                                    
        alarmed   = false;                                                   

 
        v_HL = (low_pt - high_pt);                                           
        v_HL = v_HL / max(norm(v_HL), eps);                                   

        for k = 1:steps_total                                                
 
            if trial==1 && mod(k, sample_step)==0                             
                si = si + 1;                                               
                Xs(si) = x(1);                                               
                Os(si) = x(2);                                               
            end

 
            f1 = force_MO(t, x, N);                                           
            W  = sqrt(dt) * randn(2,1);                                    
            x_pred = x + f1*dt + sigma.*W;                                   
            f2 = force_MO(t+dt, x_pred, N);                                  
            x = x + 0.5*(f1+f2)*dt + sigma.*W;                                
            t = t + dt;                                                     

 
            x(1) = max(x(1), 0);                                             
            x(2) = max(x(2), realmin);                                        

 
            if k <= burn_steps                                                
 
                pct_done = k / steps_total;                                  
                if pct_done >= next_print                                       
                    fprintf('[N=%e | trial=%d] 进度: %.0f%%（热身中）\n', ...
                        N, trial, 100*pct_done);                                
                    next_print = next_print + progress_print_every;             
                end
                continue;                                                      
            end

 
            dH = norm(x - high_pt);                                          
            dL = norm(x - low_pt );                                          
            dir_cos = dot( (x - high_pt)/max(dH,eps), v_HL );                

 
            if in_high                                                       
                if dH > r_out_high                                           
                    in_high = false;                                         
                else                                                          
                    t_high = t_high + dt;                                      
                end
            else                                                              
                if dH <= r_in_high                                            
                    in_high = true;                                            
                    t_high = t_high + dt;                                     
                end
            end

 
            if (~alarmed) && (dH >= r_cross) && (dir_cos >= cone_cos_min)     
                t_alarm = t;                                                  
                alarmed = true;                                               
            end

 
            if use_buffer_band                                               
                if abs(dH - dL) < (buf_w/2)                                  
                    t_buffer = t_buffer + dt;                                 
                else                                                        
                    if dH < dL                                                
                        t_near_H = t_near_H + dt;                           
                    else                                                       
                        t_near_L = t_near_L + dt;                              
                    end
                end
            end

 
            pct_done = k / steps_total;                                      
            if pct_done >= next_print                                         
                elapsed_steps = k - burn_steps;                               
                elapsed_time  = elapsed_steps * dt;                          
                online_frac_H = 100 * (t_high / max(elapsed_time, eps));       
                if ~alarmed                                                   
                    fprintf('[N=%e | trial=%d] 进度: %.0f%% | 在线驻留(高)=%.2f%% | 未触发\n', ...
                        N, trial, 100*pct_done, online_frac_H);                
                else                                                          
                    fprintf('[N=%e | trial=%d] 进度: %.0f%% | 在线驻留(高)=%.2f%% | 触发@t=%.3e\n', ...
                        N, trial, 100*pct_done, online_frac_H, t_alarm);       
                end
                next_print = next_print + progress_print_every;                
            end

 
 
        end

 
        analyzed_time = (steps_total - burn_steps) * dt;                        
        high_stay_trials(trial)  = t_high;                                     
        near_low_trials(trial)   = t_near_L;                                  
        near_high_trials(trial)  = t_near_H;                                  
        buffer_trials(trial)     = t_buffer;                                     
        total_time_trials(trial) = analyzed_time;                             
        if alarmed                                                              
            FPT_trials(trial) = t_alarm;                                         
            hit_trials(trial) = true;                                            
        end

 
        if trial==1                                                              
            Xs = Xs(1:si); Os = Os(1:si);                                      
            figure('Name',sprintf('Trajectory (High->Low) N=%.2e',N));         
            hold on; grid on; box on;                                           
            plot(high_pt(1), high_pt(2), 'kp', 'MarkerSize',10,'MarkerFaceColor','k','DisplayName','高稳点'); 
            plot(sad_pt(1),  sad_pt(2),  'ko', 'MarkerSize',8, 'MarkerFaceColor','y','DisplayName','鞍点');   
            plot(low_pt(1),  low_pt(2),  'ks', 'MarkerSize',8, 'MarkerFaceColor','w','DisplayName','低稳点'); 
            th = linspace(0,2*pi,360);                                        
            plot(high_pt(1)+r_in_high*cos(th),  high_pt(2)+r_in_high*sin(th),  'r--','LineWidth',1.5,'DisplayName','r_{in}');  
            plot(high_pt(1)+r_out_high*cos(th), high_pt(2)+r_out_high*sin(th), 'm--','LineWidth',1.5,'DisplayName','r_{out}'); 
            plot(high_pt(1)+r_cross*cos(th),    high_pt(2)+r_cross*sin(th),    'b--','LineWidth',1.5,'DisplayName','越界 r_{cross}'); 
            plot(Xs, Os, '-', 'LineWidth',1.2, 'DisplayName','稀疏轨迹');      
            xlabel('M'); ylabel('O');                                          
            title(sprintf('相图：从高态出发（N=%.2e）', N));                       
            legend('Location','best');                                          
            hold off;                                                         
        end
    end

 
    mean_total_time       = mean(total_time_trials);                            
    frac_high_all(i)      = 100 * mean(high_stay_trials) / mean_total_time;     
    frac_low_near_all(i)  = 100 * mean(near_low_trials)  / mean_total_time;    
    frac_high_near_all(i) = 100 * mean(near_high_trials) / mean_total_time;     
    frac_buffer_all(i)    = 100 * mean(buffer_trials)    / mean_total_time;    

 
    hit_flag_all(i) = any(hit_trials);                                      
    if any(hit_trials)                                                        
        FPT_down_all(i) = median(FPT_trials(hit_trials));                   
    else                                                                        
        FPT_down_all(i) = NaN;                                                
    end

 
    fprintf('N=%e | 命中率=%.1f%% | FPT_down(中位)=%.3e | 驻留(高)=%.1f%% | 近高=%.1f%% | 近低=%.1f%% | 缓冲=%.1f%%\n', ...
        N, 100*mean(hit_trials), FPT_down_all(i), frac_high_all(i), ...
        frac_high_near_all(i), frac_low_near_all(i), frac_buffer_all(i));       
end

 

xx = log10(N_values);                                                                    

 
figure('Name','FPT_{down} vs N（含临界线）');                                           
hold on; grid on;                                                                        
plot(xx, FPT_down_all, 'o-','LineWidth',2,'DisplayName','MFPT');          
xline(log10(crit1), 'r--', 'LineWidth',1.5, 'DisplayName','LCP');                    
xline(log10(crit2), 'g--', 'LineWidth',1.5, 'DisplayName','HCP');                   
xlabel('log_{10}(N)');                                                                
ylabel('MFPT');                                                                          
title('MFPT (High → Low State) vs. N');                                          
legend('Location','best');                                                              
hold off;                                                                                

 
figure('Name','驻留(高) vs N（含临界线）');                                                
hold on; grid on;                                                                         
plot(xx, frac_high_all, 'o-','LineWidth',2,  'DisplayName','Duration time(High)');                    
 
 
 
 
xline(log10(crit1), 'r--', 'LineWidth',1.5, 'DisplayName','LCP');                    
xline(log10(crit2), 'g--', 'LineWidth',1.5, 'DisplayName','HCP');                   
xlabel('log_{10}(N)');                                                                    
ylabel('Percentage (%)');                                                                       
title('Duration time(High → Low State) vs. N');                              
legend('Location','best');                                                                
hold off;                                                                                  

 
function f = force_MO(~, x, N)                                                           
 
    r  = 3e11;                                                                            
    s  = 2.03e-5;                                                                         
    dr = 1.36e19;                                                                       
    dd = 2.73e17;                                                                         
    a1 = 0.003; a2 = -0.1655; a3 = 3.2305; a4 = -25.8343; a5 = 71.6;                 

 
    O = max(x(2), realmin);                                                                
    delta = O / (dd + O);                                                                
    gamma = O / (dr + O);                                                                 
    omega = (1 - gamma) * (1 - delta);                                                 
    p = log10(O);                                                                       
    fai = 10^(a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5);                                      

 
    M = max(x(1), 0);                                                                       
    f1 = 0.5*omega*(N + r) - s*M - 0.5*fai*M^0.7;                                          
    f2 = omega*N - (1-omega)*r - s*M - fai*M^0.7;                                        
    f  = [f1; f2];                                                                         
end                                                                                      
