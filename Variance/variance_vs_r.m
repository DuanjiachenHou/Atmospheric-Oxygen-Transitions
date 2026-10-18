clear;                       
clc;                         

stdata = [  
    1*10^9.6, 1e15,    NaN,         NaN, 3.60592e14, 1.44193e19;  
    1*10^9.7, 1e15,    NaN,         NaN, 3.90566e14, 1.39204e19;  
    1*10^9.8, 1e15,    NaN,         NaN, 4.30739e14, 1.33217e19;  
    1*10^9.9, 1e15,    NaN,         NaN, 4.85307e14, 1.26148e19;  
    1*10^10.0,1e15,    NaN,         NaN, 5.60103e14, 1.17988e19;  
    1*10^10.1,1e15,    NaN,         NaN, 6.62756e14, 1.08836e19;  
    1*10^10.2,1e15,    NaN,         NaN, 8.02910e14, 9.89514e18;  
    1*10^10.3,1e15,    NaN,         NaN, 9.91630e14, 8.87484e18;  
    1*10^10.4,1e15,    NaN,         NaN, 1.23966e15, 7.86910e18;  
    1*10^10.5,1e15,    NaN,         NaN, 1.55810e15, 6.91370e18;  
    1*10^10.6,1e15,    NaN,         NaN, 1.96114e15, 6.02545e18;  
    1*10^10.7, 1e15,  2.46890e15, 1.12249e14, 2.46890e15, 5.20691e18;  
    1*10^10.8, 1e15,  3.10816e15, 6.11696e13, 3.10816e15, 4.45592e18;  
    1*10^10.9, 1e15,  3.91295e15, 3.88163e13, 3.91295e15, 3.76972e18;  
    1*10^11.0, 1e15,  4.92451e15, 2.65580e13, 4.92611e15, 3.14555e18;  
    1*10^11.1, 1e15,  6.20160e15, 1.90155e13, 6.20160e15, 2.58054e18;  
    1*10^11.2, 1e15,  7.80734e15, 1.40705e13, 7.80736e15, 2.07158e18;  
    1*10^11.3, 1e15,  9.82826e15, 1.06665e13, 9.82888e15, 1.61498e18;  
    1*10^11.4, 1e15,  1.23737e16, 8.23841e12, 1.23738e16, 1.20563e18;  
    1*10^11.5, 1e15,  1.55777e16, 6.46087e12, 1.55777e16, 8.34188e17;  
    1*10^11.6, 1e15,  1.96112e16, 5.13087e12, 1.96112e16, 4.66167e17;  
    1*10^11.7, 1e15,  2.46890e16, 4.11792e12,    NaN,         NaN;  
    1*10^11.8, 1e15,  3.10816e16, 3.33490e12,    NaN,         NaN;  
    1*10^11.9, 1e15,  3.91295e16, 2.72193e12,    NaN,         NaN;  
    1*10^12.0, 1e15,  4.92611e16, 2.23681e12,    NaN,         NaN;  
    1*10^12.1, 1e15,  6.20160e16, 1.84919e12,    NaN,         NaN;  
    1*10^12.2, 1e15,  7.80736e16, 1.53688e12,    NaN,         NaN;  
    1*10^12.3, 1e15,  9.82888e16, 1.28335e12,    NaN,         NaN;  
    1*10^12.4, 1e15,  1.23738e17, 1.07616e12,    NaN,         NaN;  
    1*10^12.5, 1e15,  1.55777e17, 9.05821e11,    NaN,         NaN;  
    1*10^12.6, 1e15,  1.96112e17, 7.65007e11,    NaN,         NaN;  
    1*10^12.7, 1e15,  2.46890e17, 6.48020e11,    NaN,         NaN;  
];  

r_vals    = stdata(:,1);    
 
x_low     = stdata(:,3);   
y_low     = stdata(:,4);     
x_high    = stdata(:,5);     
y_high    = stdata(:,6);     

 
h     = 1;                   
t_end = 1e7;              
t     = 0:h:t_end;           
Nt    = length(t);           
D     = 1*10^20.4;              

 
var_M_low  = NaN(size(r_vals));   
var_O_low  = NaN(size(r_vals));   
var_M_high = NaN(size(r_vals)); 
var_O_high = NaN(size(r_vals));   

 
has_low  = ~isnan(x_low)  & ~isnan(y_low);  
has_high = ~isnan(x_high) & ~isnan(y_high);  

total = length(r_vals);   
 
for k = 1:total  
    para = r_vals(k);     

 
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

        fprintf('进度 %2d/%2d: r=%.2e 低稳态  Var(M)=%.3e, Var(O)=%.3e\n', ...  
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

        fprintf('进度 %2d/%2d: r=%.2e 高稳态  Var(M)=%.3e, Var(O)=%.3e\n', ...  
                k, total, para, var_M_high(k), var_O_high(k));  
    end  
end  

 
 
figure;  
semilogx(r_vals(has_high), var_O_high(has_high), '-o', 'DisplayName','Var(O)\_high');  
xlabel('还原剂输入通量 r（mol/年）');  
ylabel('氧气方差 Var(O)');  
title('高稳点的氧气方差 vs r');  
 
crit1 = 1*10^10.6127;  
crit2 = 1*10^11.6303;   

 
xline(crit1, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');     
xline(crit2, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');     

 
legend('Location','best');
grid on;  

 
figure;  
semilogx(r_vals(has_high), var_M_high(has_high), '-d','Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980], 'DisplayName','Var(M)\_high');  
xlabel('还原剂输入通量 r（mol/年）');  
ylabel('甲烷方差 Var(M)');  
title('高稳点的甲烷方差 vs r');  
 
crit1 = 1*10^10.6127;   
crit2 = 1*10^11.6303;   

 
xline(crit1, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');      
xline(crit2, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');     

 
legend('Location','best');
grid on;  

 
figure;  
semilogx(r_vals(has_low), var_O_low(has_low), '-s', 'DisplayName','Var(O)\_low');  
xlabel('还原剂输入通量 r（mol/年）');  
ylabel('氧气方差 Var(O)');  
title('低稳点的氧气方差 vs r');  
 
crit1 = 1*10^10.6127;   
crit2 = 1*10^11.6303;   
 
xline(crit1, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');      
xline(crit2, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');      

 
legend('Location','best');
grid on;  

 
figure;  
 semilogx(r_vals(has_low), var_M_low(has_low), '-^', 'Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'DisplayName','Var(M)\_low'); 
xlabel('还原剂输入通量 r（mol/年）');  
ylabel('甲烷方差 Var(M)');  
title('低稳点的甲烷方差 vs r');  
 
crit1 = 1*10^10.6127;  
crit2 = 1*10^11.6303; 

 
xline(crit1, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');     
xline(crit2, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');     

 
legend('Location','best');
grid on;  


 
crit1 = 1*10^10.6127;  
crit2 = 1*10^11.6303;  

 
idx_O = ~isnan(var_O_low);            
figure;
semilogx(r_vals(idx_O), var_O_low(idx_O), '-s', 'DisplayName','Var(O)\_low');
hold on;
 
 
xline(crit1, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');    
xline(crit2, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');    
hold off;
xlabel('还原剂输入通量 r（mol/年, 对数刻度）');
ylabel('氧气方差 Var(O)');
title('低稳点的氧气方差 vs r（已跳过 NaN）');
legend('Location','best');
grid on;

 
idx_M = ~isnan(var_M_low);            
figure;
semilogx(r_vals(idx_M), var_M_low(idx_M), '-^', ...
    'Color',[0.8500 0.3250 0.0980], ...
    'MarkerFaceColor',[0.8500 0.3250 0.0980], ...
    'DisplayName','Var(M)\_low');
hold on;
xline(crit1, 'g--', 'LineWidth',1.5, ...  
    'DisplayName','高氧临界');    
xline(crit2, 'r--', 'LineWidth',1.5, ...  
    'DisplayName','低氧临界');      
hold off;
xlabel('还原剂输入通量 r（mol/年, 对数刻度）');
ylabel('甲烷方差 Var(M)');
title('低稳点的甲烷方差 vs r（已跳过 NaN）');
legend('Location','best');
grid on;
