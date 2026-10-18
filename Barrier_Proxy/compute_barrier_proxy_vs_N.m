function compute_barrier_proxy_vs_N()                                 
 

 
EPS_POS = 1e-12;                                                                      
r  = 3e11;                                                                           
s  = 2.03e-5; dr = 1.36e19; dd = 2.73e17;                                           
a1 = 0.003; a2 = -0.1655; a3 = 3.2305; a4 = -25.8343; a5 = 71.6;                     

 
N_list = [10^14.9, 10^15.0, 10^15.1, 10^15.2, 10^15.3, 10^15.4, 10^15.5, 10^15.6];  
LOW  = [ 1.47784E+16, 1.47783E+16, 1.47781E+16, 1.47781E+16, 1.47782E+16, 1.47784E+16, 1.47784E+16, 1.47785E+16; ... 
         4.90558E+12, 6.82215E+12, 9.74415E+12, 1.44235E+13, 2.24673E+13, 3.79197E+13, 7.46981E+13, 2.86412E+14];    
HIGH = [ 1.47783E+16, 1.47783E+16, 1.47783E+16, 1.47783E+16, 1.47783E+16, 1.47783E+16, 1.47783E+16, 1.47783E+16; ... 
         3.72183E+17, 9.16684E+17, 1.46898E+18, 2.11124E+18, 2.86041E+18, 3.72684E+18, 4.71910E+18, 5.84536E+18];   
SAD  = [ 1.17877e+16, 1.44261e+16, 1.3815e+16, 1.582131e+16, 1.517638e+16, 1.449530e+16, 1.52831e+16, 1.50176e+16; 
         9.01337e+16, 5.4242e+16, 3.26425e+16, 1.822074e+16, 1.198135e+16, 4.915387e+15, 2.8257e+15, 7.33669e+14]; 

 
tol_hit_rel = 2e-3;                                                                    
ds_rel      = 0.01;                                                                     
Kmax        = 4000;                                                                     
shrink_max  = 12;                                                                      
ns_line     = 400;                                                                      

 
dU_low  = nan(size(N_list)); hit_low  = false(size(N_list));                             
dU_high = nan(size(N_list)); hit_high = false(size(N_list));                             
fprintf('—— 沿 x 空间单位步长反向流积分势垒（非负、非零、含兜底） ——\n');                  

for i = 1:numel(N_list)                                                                 
    tic;                                                                               
    N  = N_list(i);                                                                     
    xL = LOW(:,i); xH = HIGH(:,i); xS = SAD(:,i);                                      

 
    Jx = jacobian_XY_x(xS, N, r, s, dr, dd, a1,a2,a3,a4,a5);                           
    [V,D] = eig(Jx); lam = real(diag(D));                                               
    [~,idx_stab] = min(lam);                                                              
    vS = real(V(:,idx_stab)); vS = vS / max(norm(vS),eps);                                
    dirL = sign( (xL - xS)' * vS ); if dirL==0, dirL = +1; end                           
    dirH = sign( (xH - xS)' * vS ); if dirH==0, dirH = -1; end                          
    eps_x = 1e-3 * norm(xL - xS);                                                          
    if eps_x==0, eps_x = 1e-3 * norm(xH - xS); end                                         

 
    [dU_low(i), hit_low(i)]  = follow_and_accumulate(xS + dirL*eps_x*vS, xL, N, r, s, dr, dd, ...
                                    a1,a2,a3,a4,a5, ds_rel, Kmax, shrink_max, tol_hit_rel, EPS_POS);

 
    [dU_high(i), hit_high(i)] = follow_and_accumulate(xS + dirH*eps_x*vS, xH, N, r, s, dr, dd, ...
                                    a1,a2,a3,a4,a5, ds_rel, Kmax, shrink_max, tol_hit_rel, EPS_POS); 

 
    if ~hit_low(i)                                                                         
        dU_low(i)  = max(line_positive_projection(xL, xS, N, r, s, dr, dd, a1,a2,a3,a4,a5, ns_line), EPS_POS); 
        warning('N=10^{%.2f}: 低稳支未命中，用直线正投影兜底。', log10(N));                            
    end
    if ~hit_high(i)                                                                         
        dU_high(i) = max(line_positive_projection(xH, xS, N, r, s, dr, dd, a1,a2,a3,a4,a5, ns_line), EPS_POS); 
        warning('N=10^{%.2f}: 高稳支未命中，用直线正投影兜底。', log10(N));                           
    end

 
    tagL = tern(hit_low(i),'命中','兜底'); tagH = tern(hit_high(i),'命中','兜底');                  
    fprintf('N=10^{%.2f}: ΔU_low=%8.3e (%s)   ΔU_high=%8.3e (%s) | 用时 %.2fs\n', ...
        log10(N), dU_low(i), tagL, dU_high(i), tagH, toc);                                         
end

 

 
figure;                                                                               
plot(log10(N_list), dU_low,  'o-','LineWidth',2); grid on;                            

crit1 = 10^14.8999;                                                                    
crit2 = 10^15.6026;                                                                
xline(log10(crit1), 'r--', 'LineWidth',1.5, 'DisplayName','低氧临界');               
xline(log10(crit2), 'g--', 'LineWidth',1.5, 'DisplayName','高氧临界');                
xlabel('log_{10}(N)','FontSize',12);                                                
ylabel('\Delta U（沿轨迹，低稳 \rightarrow 鞍）','FontSize',12);                       
title('沿 x 空间轨迹的势垒：低稳 \rightarrow 鞍（非负、非零）','FontSize',14);          

 
figure;                                                                                 
plot(log10(N_list), dU_high, 's-','LineWidth',2, 'DisplayName','\Delta U_{high}'); hold on; 
grid on;                                                                              

crit1 = 10^14.8999;                                                                   
crit2 = 10^15.6026;                                                                   
xline(log10(crit1), 'r--', 'LineWidth',1.5, 'DisplayName','低氧临界');                
xline(log10(crit2), 'g--', 'LineWidth',1.5, 'DisplayName','高氧临界');                 

xlabel('log_{10}(N)','FontSize',12);                                                
ylabel('\Delta U（沿轨迹，高稳 \rightarrow 鞍）','FontSize',12);                      
title('沿 x 空间轨迹的势垒：高稳 \rightarrow 鞍（非负、非零）','FontSize',14);          
legend('show','Location','best');                                                     

end                                                                                    

 

function [U, hit] = follow_and_accumulate(x0, xT, N, r, s, dr, dd, a1,a2,a3,a4,a5, ...
                                          ds_rel, Kmax, shrink_max, tol_hit_rel, EPS_POS)  
    U   = 0;                                                                               
    hit = false;                                                                         
    Lref= norm(xT - x0); if Lref==0, Lref = max(norm(xT),1); end                          
    ds0 = ds_rel * Lref;                                                                  
    x   = x0;                                                                            
    for k = 1:Kmax                                                                         
        f  = force_XY_x(x, N, r, s, dr, dd, a1,a2,a3,a4,a5);                              
        nf = norm(f); if nf < 1e-300, break; end                                          
        dir = - f / nf;                                                                  
        ds  = ds0;                                                                         
        ok  = false;                                                                       
        for sh = 1:shrink_max                                                           
            x_try = x + ds * dir;                                                         
            if all(isfinite(x_try)) && x_try(2) > 0 && x_try(1) > 0                      
                ok = true; break;                                                       
            else
                ds = 0.5*ds;                                                             
            end
        end
        if ~ok, break; end                                                               
        xm = 0.5*(x + x_try);                                                           
        fm = force_XY_x(xm, N, r, s, dr, dd, a1,a2,a3,a4,a5);                             
        seg = max(0, norm(fm) * norm(x_try - x));                                      
        U   = U + seg;                                                                   
        x   = x_try;                                                                   
        if rel_log_dist(x, xT) < tol_hit_rel                                           
            hit = true;                                                                    
            break;                                                                         
        end
    end
    U = max(U, EPS_POS);                                                               
end

function d = rel_log_dist(x, a)                                                          
    lx = [log10(max(x(1),eps)); log10(max(x(2),eps))];                                   
    la = [log10(max(a(1),eps)); log10(max(a(2),eps))];                                    
    d  = norm(lx - la)/max(norm(la),1);                                                   
end

function dU = line_positive_projection(xA, xB, N, r, s, dr, dd, a1,a2,a3,a4,a5, ns)        
    e = xB - xA; L = norm(e); if L==0, dU=0; return; end                                  
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
    M = XY(1,:); O = max(XY(2,:), eps);                                                 
    delta = O./(dd + O);                                                                   
    gamma = O./(dr + O);                                                                   
    omega = (1 - gamma).*(1 - delta);                                                     
    p  = log10(O);                                                                        
    ap = a1*p.^4 + a2*p.^3 + a3*p.^2 + a4*p + a5;                                        
    fai= 10.^ap;                                                                           
    Fx = 0.5*omega.*(N + r) - s.*M - 0.5*fai.*M.^0.7;                                      
    Fy =      omega.*N       - (1 - omega).*r - s.*M -     fai.*M.^0.7;                   
    f  = [Fx; Fy];                                                                        
end

function J = jacobian_XY_x(x, N, r, s, dr, dd, a1, a2, a3, a4, a5)                      
    M = x(1); O = max(x(2),eps);                                                           
    delta = O/(dd + O);    d_delta = dd/(dd + O)^2;                                      
    gamma = O/(dr + O);    d_gamma = dr/(dr + O)^2;                                     
    omega = (1 - gamma)*(1 - delta);                                                     
    d_omega = -(1 - delta)*d_gamma - (1 - gamma)*d_delta;                                
    p  = log10(O);                                                                      
    ap = a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5;                                           
    fai= 10^ap;                                                                            
    dap_dp  = 4*a1*p^3 + 3*a2*p^2 + 2*a3*p + a4;                                         
    dfai_dO = (fai/O)*dap_dp;                                                            
    dFx_dM = -s - 0.5*fai*0.7*M^(-0.3);                                                 
    dFy_dM = -s -     fai*0.7*M^(-0.3);                                                   
    dFx_dO = 0.5*(N + r)*d_omega - 0.5*dfai_dO*M^0.7;                                  
    dFy_dO =       (N + r)*d_omega -     dfai_dO*M^0.7;                                
    J = [dFx_dM, dFx_dO; dFy_dM, dFy_dO];                                              
end

function y = tern(cond, a, b)                                                           
    if cond, y = a; else, y = b; end                                                    
end
