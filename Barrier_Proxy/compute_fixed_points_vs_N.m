 
clear; clc; close all;                                                             

 
N  = 1*10^15.0;                                                                    
r  = 3e11;                                                                         
s  = 2.03e-5; dr = 1.36e19; dd = 2.73e17;                                          
a1 = 0.003; a2 = -0.1655; a3 = 3.2305; a4 = -25.8343; a5 = 71.6;                   
dom_u = [15, 17.5];                                                                  
dom_v = [11, 19];                                                                    
nu = 120; nv = 120;                                                                 

 
odef = @(t,x) force_MO(x, N, r, s, dr, dd, a1, a2, a3, a4, a5);                   
opts_ode = odeset('RelTol',1e-6,'AbsTol',1e-9);                                     
x_init_high = [1e16;   1e18];                                                        
x_init_low  = [1.5e16; 1.6e14];                                                      
[~,Xh] = ode23s(@(t,x) odef(t,x), [0, 5e12], x_init_high, opts_ode);               
[~,Xl] = ode23s(@(t,x) odef(t,x), [0, 5e12], x_init_low,  opts_ode);                
stab_high0 = Xh(end,:)';                                                          
stab_low0  = Xl(end,:)';                                                             

 
toUV  = @(x) [log10(x(1)); log10(x(2))];                                             
toXO  = @(z) [10.^z(1); 10.^z(2)];                                                   
Fuv   = @(z) force_XY(toXO(z), N, r, s, dr, dd, a1, a2, a3, a4, a5);                 
opts  = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12, ...
                     'StepTolerance',1e-12,'MaxFunctionEvaluations',2e4);            
[zH,fh,flH] = fsolve(Fuv, toUV(stab_high0), opts);                                   
[zL,fl,flL] = fsolve(Fuv, toUV(stab_low0),  opts);                                   
stab_high = toXO(zH);                                                                
stab_low  = toXO(zL);                                                                
fprintf('精化后稳点：高稳 M=%e, O=%e；低稳 M=%e, O=%e\n', ...
        stab_high(1),stab_high(2),stab_low(1),stab_low(2));                          

 
 
u = linspace(dom_u(1), dom_u(2), nu);                                               
v = linspace(dom_v(1), dom_v(2), nv);                                               
[UU,VV] = meshgrid(u,v);                                                             
MM = 10.^UU;  OO = 10.^VV;                                                          
[Fx,Fy] = force_XY([MM(:)'; OO(:)'], N, r, s, dr, dd, a1, a2, a3, a4, a5);           
Fx = reshape(Fx, nv, nu);  Fy = reshape(Fy, nv, nu);                                 
signFx = sign(Fx); signFy = sign(Fy);                                               
cand_uv = [];                                                                        
for i = 1:nv-1                                                                     
    for j = 1:nu-1                                                                   
        hasZeroFx = (signFx(i,j)*signFx(i,j+1)<=0) || (signFx(i+1,j)*signFx(i+1,j+1)<=0); 
        hasZeroFy = (signFy(i,j)*signFy(i,j+1)<=0) || (signFy(i+1,j)*signFy(i+1,j+1)<=0); 
        if hasZeroFx && hasZeroFy                                                    
            uc = mean([u(j), u(j+1)]);                                              
            vc = mean([v(i), v(i+1)]);                                             
            cand_uv = [cand_uv, [uc; vc]];                                          
        end
    end
end
 
e  = (stab_high - stab_low); e = e / norm(e);                                       
L  = norm(stab_high - stab_low);                                                     
g  = @(s) force_MO(stab_low + s*e, N, r, s, dr, dd, a1, a2, a3, a4, a5)'*e;         
ss = linspace(0, L, 100);                                                           
gs = arrayfun(@(s0) g(s0), ss);                                                      
idx = find(gs(1:end-1).*gs(2:end) <= 0, 1, 'first');                                
if ~isempty(idx)                                                                     
    a = ss(idx); b = ss(idx+1);                                                      
    for it = 1:30                                                                    
        c = 0.5*(a+b);                                                              
        if g(a)*g(c) <= 0, b = c; else, a = c; end                                  
    end
    s_star = 0.5*(a+b);                                                              
    x_star = stab_low + s_star*e;                                                   
    cand_uv = [cand_uv, toUV(x_star)];                                              
end
 
if ~isempty(cand_uv)                                                                
    cand_uv = unique(round(cand_uv', 6, 'significant'),'rows','stable')';            
end

 
roots = [];                                                                          
for k = 1:size(cand_uv,2)                                                            
    [zk, fk, flg] = fsolve(Fuv, cand_uv(:,k), opts);                              
    if flg==1 && all(isfinite(fk))                                                  
        xk = toXO(zk);                                                             
        if all(isfinite(xk)) && all(xk>0)                                           
            roots = [roots, xk];                                                   
        end
    end
end
 
if ~isempty(roots)                                                                  
    roots = [roots, stab_high, stab_low];                                          
    roots = unique(round(log10(roots)', 6, 'significant'),'rows','stable')';       
    roots = 10.^roots;                                                              
end

 
saddle_list = [];                                                                    
stable_list = [];                                                                     
for i = 1:size(roots,2)                                                               
    x0 = roots(:,i);                                                                  
    J  = jacobian_XY(x0, N, r, s, dr, dd, a1, a2, a3, a4, a5);                        
    lam = eig(J);                                                                    
    if det(J) < 0                                                                    
        saddle_list = [saddle_list, x0];                                              
    else                                                                             
        if all(real(lam)<0)                                                           
            stable_list = [stable_list, x0];                                          
        end
    end
end

 
if isempty(saddle_list)                                                              
    warning('当前 N 下未检测到鞍点：可能单稳态或初猜不足（建议扩大网格/检查参数区间）。'); 
else                                                                                 
 
    dists = vecnorm(saddle_list - stab_low, 2, 1) + vecnorm(saddle_list - stab_high, 2, 1); 
    [~,I] = min(dists);                                                              
    saddle = saddle_list(:,I);                                                       
    fprintf('鞍点 M=%e, O=%e\n', saddle(1), saddle(2));                             
end

 
[xg, yg] = deal(logspace(dom_u(1), dom_u(2), 400), logspace(dom_v(1), dom_v(2), 400));
[X, Y] = meshgrid(xg, yg);                                                            
[Fx2, Fy2] = force_XY([X(:)'; Y(:)'], N, r, s, dr, dd, a1, a2, a3, a4, a5);          
Fx2 = reshape(Fx2, size(X)); Fy2 = reshape(Fy2, size(Y));                            
figure; hold on; set(gca,'XScale','log','YScale','log');                              
contour(X, Y, Fx2, [0 0], 'b','LineWidth',1);                                         
contour(X, Y, Fy2, [0 0], 'r','LineWidth',1);                                         
plot(stab_high(1), stab_high(2), 'go','MarkerFaceColor','g','DisplayName','高稳点');  
plot(stab_low(1),  stab_low(2),  'mo','MarkerFaceColor','m','DisplayName','低稳点'); 
if ~isempty(saddle_list), plot(saddle_list(1,:), saddle_list(2,:), 'yp', ...
     'MarkerSize',12,'MarkerFaceColor','y','DisplayName','鞍点候选'); end             
if exist('saddle','var'), plot(saddle(1), saddle(2), 'kp','MarkerSize',14, ...
     'MarkerFaceColor','k','DisplayName','最终鞍点'); end                             
xlabel('M'); ylabel('O'); title(sprintf('N = 10^{%.2f} 的平衡点与鞍点', log10(N)));   
legend('Location','best'); grid on; axis tight;                                        
ylim([10^dom_v(1), 10^dom_v(2)]); xlim([10^dom_u(1), 10^dom_u(2)]);                  

 

function dx = force_MO(x, N, r, s, dr, dd, a1, a2, a3, a4, a5)                          
    [Fx,Fy] = force_XY(x, N, r, s, dr, dd, a1, a2, a3, a4, a5);                         
    dx = [Fx; Fy];                                                                       
end

function [Fx, Fy] = force_XY(XY, N, r, s, dr, dd, a1, a2, a3, a4, a5)                    
    M = XY(1,:); O = XY(2,:);                                                          
    delta = O./(dd + O);                                                                
    gamma = O./(dr + O);                                                               
    omega = (1 - gamma).*(1 - delta);                                                   
    p  = log10(O);                                                                       
    ap = a1*p.^4 + a2*p.^3 + a3*p.^2 + a4*p + a5;                                         
    fai= 10.^ap;                                                                           
    Fx = 0.5*omega.*(N + r) - s.*M - 0.5*fai.*M.^0.7;                                     
    Fy =      omega.*N       - (1 - omega).*r - s.*M -     fai.*M.^0.7;                    
end

function J = jacobian_XY(x, N, r, s, dr, dd, a1, a2, a3, a4, a5)                          
    M = x(1); O = x(2);                                                                  
    delta = O/(dd + O);    d_delta = dd/(dd + O)^2;                                       
    gamma = O/(dr + O);    d_gamma = dr/(dr + O)^2;                                       
    omega = (1 - gamma)*(1 - delta);                                                      
    d_omega = -(1 - delta)*d_gamma - (1 - gamma)*d_delta;                                
    p  = log10(O);                                                                          
    ap = a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5;                                           
    fai= 10^ap;                                                                            
    dap_dp = 4*a1*p^3 + 3*a2*p^2 + 2*a3*p + a4;                                           
    dfai_dO = (fai/O)*dap_dp;                                                             
    dFx_dM = -s - 0.5*fai*0.7*M^(-0.3);                                                    
    dFy_dM = -s -     fai*0.7*M^(-0.3);                                                   
    dFx_dO = 0.5*(N + r)*d_omega - 0.5*dfai_dO*M^0.7;                                      
    dFy_dO =       (N + r)*d_omega -     dfai_dO*M^0.7;                                     
    J = [dFx_dM, dFx_dO; dFy_dM, dFy_dO];                                                  
end
