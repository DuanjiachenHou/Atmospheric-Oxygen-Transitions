function compute_barrier_proxy_vs_r()                                   

 
EPS_POS = 1e-12;                                                           
N = 1e15;                                                                 
s  = 2.03e-5;                                                              
dr = 1.36e19;                                                               
dd = 2.73e17;                                                               
a1 = 0.003; a2 = -0.1655; a3 = 3.2305; a4 = -25.8343; a5 = 71.6;       

 
r_list = [1*10^10.65, 1*10^10.7, 1*10^10.8, 1*10^10.9, 1*10^11.0, ...        
          1*10^11.1, 1*10^11.2, 1*10^11.3, 1*10^11.4, 1*10^11.5, 1*10^11.6]; 

LOW = [                                                                       
    2.20041E+15, 2.46890E+15, 3.10816E+15, 3.91295E+15, 4.92451E+15, 6.20160E+15, 7.80734E+15, ...
    9.82826E+15, 1.23737E+16, 1.55777E+16, 1.96112E+16;
    1.76933E+14, 1.12249E+14, 6.11696E+13, 3.88163E+13, 2.65580E+13, 1.90155E+13, 1.40705E+13, ...
    1.06665E+13, 8.23841E+12, 6.46087E+12, 5.13087E+12
];

HIGH = [                                                                      
    2.20041E+15, 2.46890E+15, 3.10816E+15, 3.91295E+15, 4.92611E+15, 6.20160E+15, 7.80736E+15, ...
    9.82888E+15, 1.23738E+16, 1.55777E+16, 1.96112E+16;
    5.60757E+18, 5.20691E+18, 4.45592E+18, 3.76972E+18, 3.14555E+18, 2.58054E+18, 2.07158E+18, ...
    1.61498E+18, 1.20563E+18, 8.34188E+17, 4.66167E+17
];

SAD = [                                                                      
    2.18675E+15, 2.48163E+15, 3.12599E+15, 3.99488E+15, 4.89892E+15, 5.78229E+15, 7.75755E+15, ...
    9.91381E+15, 1.21328E+16, 1.55052E+16, 1.95311E+16;
    1.00000E+15, 1.78085E+15, 3.39882E+15, 5.39299E+15, 8.73373E+15, 1.19015E+16, 2.15443E+16, ...
    3.41849E+16, 4.50957E+16, 7.49353E+16, 1.56852E+17
];

 
tol_hit_log = 2e-3;                                                        
tol_hit_lin = 2e-2;                                                        
ds_rel      = 0.01;                                                         
Kmax        = 5000;                                                        
shrink_max  = 12;                                                        
ns_line     = 600;                                                  

 
dU_low  = nan(size(r_list));                                             
dU_high = nan(size(r_list));                                               
hit_low  = false(size(r_list));                                            
hit_high = false(size(r_list));                                            

fprintf('—— 沿 x 空间单位步长反向流积分势垒（优先命中，不命中则兜底） ——\n'); 

for i = 1:numel(r_list)                                                     
    tStart = tic;                                                         
    r  = r_list(i);                                                     
    xL = LOW(:,i);                                                       
    xH = HIGH(:,i);                                                      
    xS = SAD(:,i);                                                          

    Jx = jacobian_XY_x(xS, N, r, s, dr, dd, a1,a2,a3,a4,a5);               
    [V,Dm] = eig(Jx);                                                      
    lam = real(diag(Dm));                                                  
    [~,idx_stab] = min(lam);                                            
    vS = real(V(:,idx_stab));                                              
    vS = vS / max(norm(vS),eps);                                           

    eps_x = 1e-3 * max(norm(xL - xS), norm(xH - xS));                        
    eps_x = max(eps_x, 1e-6*max(norm(xS),1));                                

    x0p = xS + eps_x*vS;                                                  
    x0m = xS - eps_x*vS;                                                   

 
    [U_Lp, hit_Lp] = follow_and_accumulate(x0p, xL, N, r, s, dr, dd, ...      
        a1,a2,a3,a4,a5, ds_rel, Kmax, shrink_max, tol_hit_log, tol_hit_lin, EPS_POS); 
    [U_Lm, hit_Lm] = follow_and_accumulate(x0m, xL, N, r, s, dr, dd, ...     
        a1,a2,a3,a4,a5, ds_rel, Kmax, shrink_max, tol_hit_log, tol_hit_lin, EPS_POS); 

    if hit_Lp && hit_Lm                                                      
        dU_low(i)  = min(U_Lp, U_Lm);                                       
        hit_low(i) = true;                                                
    elseif hit_Lp                                                          
        dU_low(i)  = U_Lp;                                                  
        hit_low(i) = true;                                                  
    elseif hit_Lm                                                         
        dU_low(i)  = U_Lm;                                                
        hit_low(i) = true;                                                  
    else                                                                    
        hit_low(i) = false;                                      
    end                                                                    

 
    [U_Hp, hit_Hp] = follow_and_accumulate(x0p, xH, N, r, s, dr, dd, ...      
        a1,a2,a3,a4,a5, ds_rel, Kmax, shrink_max, tol_hit_log, tol_hit_lin, EPS_POS); 
    [U_Hm, hit_Hm] = follow_and_accumulate(x0m, xH, N, r, s, dr, dd, ...      
        a1,a2,a3,a4,a5, ds_rel, Kmax, shrink_max, tol_hit_log, tol_hit_lin, EPS_POS); 

    if hit_Hp && hit_Hm                                                      
        dU_high(i)  = min(U_Hp, U_Hm);                                      
        hit_high(i) = true;                                                 
    elseif hit_Hp                                                            
        dU_high(i)  = U_Hp;                                                  
        hit_high(i) = true;                                                  
    elseif hit_Hm                                                            
        dU_high(i)  = U_Hm;                                                 
        hit_high(i) = true;                                                  
    else                                                                    
        hit_high(i) = false;                                                
    end                                                                      

 
    if ~hit_low(i)                                                          
        dU_low(i) = max(line_positive_projection(xL, xS, N, r, s, dr, dd, ... 
            a1,a2,a3,a4,a5, ns_line), EPS_POS);                               
        warning('r=10^{%.2f}: 低稳支未命中，用直线正投影兜底。', log10(r));     
    end                                                                      

    if ~hit_high(i)                                                         
        dU_high(i) = max(line_positive_projection(xH, xS, N, r, s, dr, dd, ...
            a1,a2,a3,a4,a5, ns_line), EPS_POS);                               
        warning('r=10^{%.2f}: 高稳支未命中，用直线正投影兜底。', log10(r));      
    end                                                                      

    tagL = tern(hit_low(i),'命中','兜底');                                   
    tagH = tern(hit_high(i),'命中','兜底');                                  
    fprintf('r=10^{%.2f}: ΔU_low=%8.3e (%s)  ΔU_high=%8.3e (%s) | 用时 %.2fs\n', ... 
        log10(r), dU_low(i), tagL, dU_high(i), tagH, toc(tStart));            
end                                                                          

 
fprintf('\n===== 命中统计 =====\n');                                          
fprintf('low  命中率 = %.1f%%\n', 100*mean(hit_low));                        
fprintf('high 命中率 = %.1f%%\n', 100*mean(hit_high));                        
T = table(log10(r_list(:)), hit_low(:), hit_high(:), dU_low(:), dU_high(:), ... 
    'VariableNames', {'log10r','hit_low','hit_high','dU_low','dU_high'});    
disp(T);                                                                      

 
x_left  = 10.61;                                                              
x_right = 11.63;                                                              

figure;                                                                       
plot(log10(r_list), dU_low, 'o-','LineWidth',2); grid on;                     
xlabel('log_{10}(r)','FontSize',12);                                          
ylabel('\Delta U（沿轨迹，低稳 \rightarrow 鞍）','FontSize',12);               
title('沿 x 空间轨迹的势垒：低稳 \rightarrow 鞍','FontSize',14);               
xline(x_left,  'g--','LineWidth',2);                                       
xline(x_right, 'r--','LineWidth',2);                                         

figure;                                                                     
plot(log10(r_list), dU_high, 's-','LineWidth',2, 'DisplayName','\Delta U_{high}'); hold on; 
grid on;                                                                      
xlabel('log_{10}(r)','FontSize',12);                                         
ylabel('\Delta U（沿轨迹，高稳 \rightarrow 鞍）','FontSize',12);               
title('沿 x 空间轨迹的势垒：高稳 \rightarrow 鞍','FontSize',14);              
legend('show','Location','best');                                           
xline(x_left,  'g--','LineWidth',2);                                         
xline(x_right, 'r--','LineWidth',2);                                       

end                                                                          

 

function [U, hit] = follow_and_accumulate(x0, xT, N, r, s, dr, dd, ...        
    a1,a2,a3,a4,a5, ds_rel, Kmax, shrink_max, tol_hit_log, tol_hit_lin, EPS_POS) 

U   = 0;                                                                  
hit = false;                                                                
Lref = norm(xT - x0);                                                       
if Lref==0, Lref = max(norm(xT),1); end                                     
ds0  = ds_rel * Lref;                                                        
x    = x0;                                                                   

for k = 1:Kmax                                                              
    f  = force_XY_x(x, N, r, s, dr, dd, a1,a2,a3,a4,a5);                     
    nf = norm(f);                                                          
    if nf < 1e-300, break; end                                          
    dir = -f / nf;                                                         

    ds  = ds0;                                                               
    ok  = false;                                                          

    for sh = 1:shrink_max                                                  
        x_try = x + ds * dir;                                               
        if all(isfinite(x_try)) && x_try(1) > 0 && x_try(2) > 0             
            ok = true;                                                       
            break;                                                        
        else
            ds = 0.5*ds;                                                    
        end
    end

    if ~ok, break; end                                                       

    xm = 0.5*(x + x_try);                                                   
    fm = force_XY_x(xm, N, r, s, dr, dd, a1,a2,a3,a4,a5);                    
    dx = x_try - x;                                                          

    seg = max(0, -(fm(1)*dx(1) + fm(2)*dx(2)));                             
    U   = U + seg;                                                          
    x   = x_try;                                                           

    if is_hit(x, xT, tol_hit_log, tol_hit_lin)                                
        hit = true;                                                         
        break;                                                              
    end
end

U = max(U, EPS_POS);                                                        

end                                                                         

function tf = is_hit(x, a, tol_log, tol_lin)                               
dlog = rel_log_dist(x, a);                                                   
dlin = rel_lin_dist(x, a);                                                  
tf = (dlog < tol_log) || (dlin < tol_lin);                                   
end                                                                          

function d = rel_log_dist(x, a)                                               
lx = [log10(max(x(1),eps)); log10(max(x(2),eps))];                            
la = [log10(max(a(1),eps)); log10(max(a(2),eps))];                           
d  = norm(lx - la) / max(norm(la),1);                                  
end                                                                         

function d = rel_lin_dist(x, a)                                               
d = norm(x - a) / max(norm(a),1);                                            
end                                                                          

function dU = line_positive_projection(xA, xB, N, r, s, dr, dd, ...            
    a1,a2,a3,a4,a5, ns)

e = xB - xA;                                                                
L = norm(e);                                                                  
if L==0, dU=0; return; end                                                   
e = e / L;                                                                   

if mod(ns,2)==1, ns = ns+1; end                                             
t  = linspace(0,1,ns+1);                                                     
XY = xA + e*(L*t);                                                       
F  = force_XY_x(XY, N, r, s, dr, dd, a1,a2,a3,a4,a5);                        
val = max(0, -(F(1,:).*e(1) + F(2,:).*e(2)));                                

h = 1/(numel(t)-1);                                                          
S = val(1) + val(end) + 4*sum(val(2:2:end-1)) + 2*sum(val(3:2:end-2));      
dU = L * (h/3) * S;                                                          

end                                                                           

function f = force_XY_x(XY, N, r, s, dr, dd, a1, a2, a3, a4, a5)            
M = XY(1,:);                                                                  
O = max(XY(2,:), eps);                                                       

delta = O./(dd + O);                                                          
gamma = O./(dr + O);                                                          
omega = (1 - gamma).*(1 - delta);                                           

p  = log10(O);                                                             
ap = a1*p.^4 + a2*p.^3 + a3*p.^2 + a4*p + a5;                                
fai = 10.^ap;                                                               

Fx = 0.5*omega.*(N + r) - s.*M - 0.5*fai.*M.^0.7;                             
Fy =      omega.*N       - (1 - omega).*r - s.*M -     fai.*M.^0.7;          

f  = [Fx; Fy];                                                                
end                                                                         

function J = jacobian_XY_x(x, N, r, s, dr, dd, a1, a2, a3, a4, a5)            
M = x(1);                                                                    
O = max(x(2),eps);                                                            

delta   = O/(dd + O);                                                      
d_delta = dd/(dd + O)^2;                                                     

gamma   = O/(dr + O);                                                      
d_gamma = dr/(dr + O)^2;                                                    

omega   = (1 - gamma)*(1 - delta);                                          
d_omega = -(1 - delta)*d_gamma - (1 - gamma)*d_delta;                         

p  = log10(O);                                                               
ap = a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5;                                    
fai = 10^ap;                                                                  

dap_dp  = 4*a1*p^3 + 3*a2*p^2 + 2*a3*p + a4;                                  
dfai_dO = (fai/O)*dap_dp;                                                     

dFx_dM = -s - 0.5*fai*0.7*M^(-0.3);                                           
dFy_dM = -s -     fai*0.7*M^(-0.3);                                           

dFx_dO = 0.5*(N + r)*d_omega - 0.5*dfai_dO*M^0.7;                            
dFy_dO =      (N + r)*d_omega -     dfai_dO*M^0.7;                             

J = [dFx_dM, dFx_dO; dFy_dM, dFy_dO];                                        
end                                                                           

function y = tern(cond, a, b)                                               
if cond                                                                       
    y = a;                                                                    
else                                                                          
    y = b;                                                                   
end                                                                         
end                                                                           
