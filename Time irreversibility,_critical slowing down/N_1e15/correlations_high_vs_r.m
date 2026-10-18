function demo_highbranch_from_table(r)                                        
 
    clear; clc; close all;                                                   

 
    data = [...
        10^(9.6), 3.60592e14, 1.44193e19;
        10^(9.7), 3.90566e14, 1.39204e19;
        10^(9.8), 4.30739e14, 1.33217e19;
        10^(9.9), 4.85307e14, 1.26148e19;
        10^(10.0),5.60103e14, 1.17988e19;
        10^(10.1),6.62756e14, 1.08836e19;
        10^(10.2),8.02910e14, 9.89514e18;
        10^(10.3),9.91630e14, 8.87484e18;
        10^(10.4),1.23966e15, 7.86910e18;
        10^(10.5),1.55810e15, 6.91370e18;
        10^(10.6),1.96114e15, 6.02545e18;
        10^(10.65),2.20041e15,5.60757e18;
        10^(10.7),2.46890e15, 5.20691e18;
        10^(10.8),3.10816e15, 4.45592e18;
        10^(10.9),3.91295e15, 3.76972e18;
        10^(11.0),4.92611e15, 3.14555e18;
        10^(11.1),6.20160e15, 2.58054e18;
        10^(11.2),7.80736e15, 2.07158e18;
        10^(11.3),9.82888e15, 1.61498e18;
        10^(11.4),1.23738e16,1.20563e18;
        10^(11.5),1.55777e16,8.34188e17;
        10^(11.6),1.96112e16,4.66167e17 ];

 
    r_list = data(:,1);  
    x_list = data(:,2);   
    y_list = data(:,3);   

 
    N_const = 1e15;       
    s       = 2.03e-5;    
    dr      = 1.36e19;   
    dd      = 2.73e17;    
    a1 = 0.003; a2 = -0.1655; a3 = 3.2305; a4 = -25.8343; a5 = 71.6;  

 
    D = [1e10, 0;
         0,    1e30];

 
    nData    = numel(r_list);           
    tau_vals = zeros(nData,1);          
    dcc_vals = zeros(nData,1);          

    for i = 1:nData
        r_curr = r_list(i);             
        M_ss   = x_list(i);              
        O_ss   = y_list(i);             

 
        A = calcJacobian(M_ss, O_ss, N_const, r_curr, s, dr, dd, a1, a2, a3, a4, a5);

 
        Sigma = lyap(A, 2*D);

 
        tau_val = calcAutoCorrTime(A);

 
        if isnan(tau_val)
            dcc_val = NaN;             
        else
            tMax   = 5 * tau_val;       
            nPts   = 300;               
            dcc_val = calcDeltaCC(A, Sigma, tMax, nPts);
        end

        tau_vals(i) = tau_val;         
        dcc_vals(i) = dcc_val;

        fprintf('i=%d: r=%.3g, tau=%.4g, ΔCC=%.4g\n', i, r_curr, tau_val, dcc_val);
    end

 
    crit1 = 1*10^10.6127;           
    crit2 = 1*10^11.6303;            
 
 
    figure;                                        
    semilogx(r_list, tau_vals, 'o-','LineWidth',1.2,'MarkerSize',6,'DisplayName','\tau');  
    hold on;                                        
    xline(crit1, 'g--','LineWidth',1.5,'DisplayName','高氧临界');     
    xline(crit2, 'r--','LineWidth',1.5,'DisplayName','低氧临界');     
    xlabel('r'); ylabel('\tau');                     
    title('高稳态：自相关时间 \tau vs r');          
    legend('Location','best');                     
    grid on;                                     

 
    figure;                                       
    semilogx(r_list, dcc_vals, 's-','LineWidth',1.2,'MarkerSize',6,'DisplayName','\Delta CC');  
    hold on;                                         
    xline(crit1, 'g--','LineWidth',1.5,'DisplayName','高氧临界');    
    xline(crit2, 'r--','LineWidth',1.5,'DisplayName','低氧临界');     
    xlabel('r'); ylabel('\Delta CC');                 
    title('高稳态：\Delta CC vs r');                
    legend('Location','best');                    
    grid on;                                  
end

 
function A = calcJacobian(M_ss, O_ss, N, r, s, dr, dd, a1, a2, a3, a4, a5)
 
    F1 = @(M,O)(0.5*omega(O,dr,dd)*(N + r) - s*M - 0.5*phi(O,a1,a2,a3,a4,a5)*M^0.7);
    F2 = @(M,O)(omega(O,dr,dd)*N - (1-omega(O,dr,dd))*r - s*M - phi(O,a1,a2,a3,a4,a5)*M^0.7);
 
    eps_rel = 1e-6; dM = max(M_ss*eps_rel,1e-8); dO = max(O_ss*eps_rel,1e-8);
 
    dF1_dM = (F1(M_ss+dM,O_ss) - F1(M_ss-dM,O_ss)) / (2*dM);
    dF1_dO = (F1(M_ss,O_ss+dO) - F1(M_ss,O_ss-dO)) / (2*dO);
    dF2_dM = (F2(M_ss+dM,O_ss) - F2(M_ss-dM,O_ss)) / (2*dM);
    dF2_dO = (F2(M_ss,O_ss+dO) - F2(M_ss,O_ss-dO)) / (2*dO);
    A = [dF1_dM, dF1_dO; dF2_dM, dF2_dO];    
end

function val = omega(O_val, dr, dd)
    gamma = O_val/(dr + O_val); delta = O_val/(dd + O_val);
    val = (1 - gamma)*(1 - delta);         
end

function val = phi(O_val, a1,a2,a3,a4,a5)
    p = log10(O_val);
    val = 10^(a1*p^4 + a2*p^3 + a3*p^2 + a4*p + a5);  
end

function tau_val = calcAutoCorrTime(A)
    reVals = real(eig(A)); neg = reVals(reVals<0);
    if isempty(neg)
        tau_val = NaN;    
    else
        tau_val = -1/max(neg);  
    end
end

function dcc_val = calcDeltaCC(A, Sigma, tMax, nPts)
    t = linspace(0,tMax,nPts); diffSq = zeros(1,nPts);
    for k = 1:nPts
        Msig = expm(A*t(k))*Sigma;
        diffSq(k) = (Msig(1,2) - Msig(2,1))^2; 
    end
    area = trapz(t,diffSq); dcc_val = sqrt(area/tMax);  
end
