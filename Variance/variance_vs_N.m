clear;                       
clc;                      

 
 
stdata = [  
     3e11, 1*10^14.0, 1.47783e16, 4.66623e11,    NaN,         NaN;
    3e11, 1*10^14.1, 1.47783e16, 5.86664e11,    NaN,         NaN;
    3e11, 1*10^14.2, 1.47783e16, 7.41629e11,    NaN,         NaN;
    3e11, 1*10^14.3, 1.47783e16, 9.43449e11,    NaN,         NaN;
    3e11, 1*10^14.4, 1.47783e16, 1.20895e12,    NaN,         NaN;
    3e11, 1*10^14.5, 1.47783e16, 1.56228e12,    NaN,         NaN;
    3e11, 1*10^14.6, 1.47783e16, 2.03885e12,    NaN,         NaN;
    3e11, 1*10^14.7, 1.47783e16, 2.69193e12,    NaN,         NaN;
    3e11, 1*10^14.8, 1.47811e16, 3.60494e12,    NaN,         NaN;
    3e11, 1*10^14.9, 1.47784e16, 4.90558e12, 1.47783e16, 3.72183e17;
    3e11, 1*10^15.0, 1.47783e16, 6.82215e12, 1.47783e16, 9.16684e17;
    3e11, 1*10^15.1, 1.47781e16, 9.74415e12, 1.47783e16, 1.46898e18;
    3e11, 1*10^15.2, 1.47781e16, 1.44235e13, 1.47783e16, 2.11124e18;
    3e11, 1*10^15.3, 1.47782e16, 2.24673e13, 1.47783e16, 2.86041e18;
    3e11, 1*10^15.4, 1.47784e16, 3.79197e13, 1.47783e16, 3.72684e18;
    3e11, 1*10^15.5, 1.47784e16, 7.46981e13, 1.47783e16, 4.71910e18;
    3e11, 1*10^15.6, 1.47785e16, 2.86412e14, 1.47783e16, 5.84536e18;
    3e11, 1*10^15.7,      NaN,         NaN, 1.47783e16, 7.11396e18;
    3e11, 1*10^15.8,      NaN,         NaN, 1.47783e16, 8.53375e18;
    3e11, 1*10^15.9,      NaN,         NaN, 1.47783e16, 1.01143e19;
    3e11, 1*10^16.0,      NaN,         NaN, 1.47783e16, 1.18657e19;
    3e11, 1*10^16.1,      NaN,         NaN, 1.47783e16, 1.37994e19;
    3e11, 1*10^16.2,      NaN,         NaN, 1.47783e16, 1.59272e19;
    3e11, 1*10^16.3,      NaN,         NaN, 1.47783e16, 1.82621e19;
    3e11, 1*10^16.4,      NaN,         NaN, 1.47783e16, 2.08180e19;
    3e11, 1*10^16.5,      NaN,         NaN, 1.47783e16, 2.36098e19;
];  

 
r_vals    = stdata(:,1);      
N_vals    = stdata(:,2);      
x_low     = stdata(:,3);     
y_low     = stdata(:,4);     
x_high    = stdata(:,5);    
y_high    = stdata(:,6);    

 
h     = 100;  %10               
t_end = 1e8;                  
t     = 0:h:t_end;          
Nt    = length(t);            
D     = 1*10^20.4;              

 
var_M_low  = NaN(size(N_vals));   
var_O_low  = NaN(size(N_vals));    
var_M_high = NaN(size(N_vals));   
var_O_high = NaN(size(N_vals));   

 
has_low  = ~isnan(x_low)  & ~isnan(y_low);  
has_high = ~isnan(x_high) & ~isnan(y_high);  

total = length(N_vals);     
 
for k = 1:total  
    para = N_vals(k);      

 
    if has_low(k)  
        x = zeros(1, Nt);   
        y = zeros(1, Nt);  
        x(1) = x_low(k);   
        y(1) = y_low(k);    

        for i = 1:Nt-1  
            [Fx, Fy] = force_XY(x(i), y(i), para);              
            x(i+1) = x(i) + Fx*h + sqrt(2*D*h)*randn();         
            y(i+1) = y(i) + Fy*h + sqrt(2*D*h)*randn();          
            if x(i+1)<0, x(i+1) = -x(i+1); end                   
            if y(i+1)<0, y(i+1) = -y(i+1); end  
        end  

        var_M_low(k) = var(x);  
        var_O_low(k) = var(y);  

        fprintf('进度 %2d/%2d: N=%.2e 低稳态  Var(M)=%.3e, Var(O)=%.3e\n', ...  
                k, total, para, var_M_low(k), var_O_low(k));  
    end  

 
    if has_high(k)  
        x = zeros(1, Nt);   
        y = zeros(1, Nt);   
        x(1) = x_high(k);   
        y(1) = y_high(k);   

        for i = 1:Nt-1  
            [Fx, Fy] = force_XY(x(i), y(i), para);  
            x(i+1) = x(i) + Fx*h + sqrt(2*D*h)*randn();  
            y(i+1) = y(i) + Fy*h + sqrt(2*D*h)*randn();  
            if x(i+1)<0, x(i+1) = -x(i+1); end  
            if y(i+1)<0, y(i+1) = -y(i+1); end  
        end  

        var_M_high(k) = var(x);  
        var_O_high(k) = var(y);  

        fprintf('进度 %2d/%2d: N=%.2e 高稳态  Var(M)=%.3e, Var(O)=%.3e\n', ...  
                k, total, para, var_M_high(k), var_O_high(k));  
    end  
end  



 
crit1 = 1*10^14.8999;            
 
crit2 = 1*10^15.6026;            

 
figure;                        
semilogx(N_vals(has_high), var_O_high(has_high), '-o', ...  
    'DisplayName','Var(O)\_high'); hold on;  
xline(crit1, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');     
xline(crit2, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');     
hold off;  
xlabel('净初级生产率 N（mol/年）');  
ylabel('氧气方差 Var(O)');  
title('高稳点的氧气方差 vs N');  
legend('Location','best');  
grid on;  

 
figure;  
semilogx(N_vals(has_high), var_M_high(has_high), '-d', 'Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'DisplayName','Var(M)\_high'); hold on;  
xline(crit1, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');      
xline(crit2, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');      
hold off;  
xlabel('净初级生产率 N（mol/年）');  
ylabel('甲烷方差 Var(M)');  
title('高稳点的甲烷方差 vs N');  
legend('Location','best');  
grid on;  

 
figure;  
semilogx(N_vals(has_low), var_O_low(has_low), '-s', 'DisplayName','Var(O)\_low'); hold on;  
xline(crit1, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');      
xline(crit2, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');      
hold off;  
xlabel('净初级生产率 N（mol/年）');  
ylabel('氧气方差 Var(O)');  
title('低稳点的氧气方差 vs N');  
legend('Location','best');  
grid on;  

 
figure;  
semilogx(N_vals(has_low), var_M_low(has_low), '-^', 'Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'DisplayName','Var(M)\_low'); hold on;  
xline(crit1, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');     
xline(crit2, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');      
hold off;  
xlabel('净初级生产率 N（mol/年）');  
ylabel('甲烷方差 Var(M)');  
title('低稳点的甲烷方差 vs N');  
legend('Location','best');  
grid on;  
