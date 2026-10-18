function demo_lowbranch_from_table(N)
 
    clear; clc; close all;                                                   

 
 
    data = [...                                                                
       1e14,       1.47783e16, 4.66623e11;
       1.25893e14, 1.47783e16, 5.86664e11;  % 10^(14.1)=1.25893e14
       1.58489e14, 1.47783e16, 7.41629e11;  % 10^(14.2)=1.58489e14
       1.99526e14, 1.47783e16, 9.43449e11;  % 10^(14.3)=1.99526e14
       2.51189e14, 1.47783e16, 1.20895e12;  % 10^(14.4)
       3.16228e14, 1.47783e16, 1.56228e12;  % 10^(14.5)
       3.98107e14, 1.47783e16, 2.03885e12;  % 10^(14.6)
       5.01187e14, 1.47783e16, 2.69193e12;  % 10^(14.7)
       6.30957e14, 1.47811e16, 3.60494e12;  % 10^(14.8)
       7.94328e14, 1.47784e16, 4.90558e12;  % 10^(14.9)
       1e15,       1.47783e16, 6.82215e12;  % 10^(15.0)
       1.25893e15, 1.47781e16, 9.74415e12;  % 10^(15.1)
       1.58489e15, 1.47781e16, 1.44235e13;  % 10^(15.2)
       1.99526e15, 1.47782e16, 2.24673e13;  % 10^(15.3)
       2.51189e15, 1.47784e16, 3.79197e13;  % 10^(15.4)
       3.16228e15, 1.47784e16, 7.46981e13;  % 10^(15.5)
       3.98107e15, 1.47785e16, 2.86412e14]; % 10^(15.6)

 
    N_list = data(:,1);                                                      
    M_list = data(:,2);                                                      
    O_list = data(:,3);                                                      

 
    r  = 3e11;  s  = 2.03e-5;  dr = 1.36e19;  dd = 2.73e17;                   
    a1 = 0.003; a2 = -0.1655; a3 = 3.2305; a4 = -25.8343; a5 = 71.6;           
    %D  = [1e24, 0; 0, 1e32];                                                  
     D  = [1e10, 0; 0, 1e30];   
 
    nData    = size(data,1);                                                  
    tau_vals = zeros(nData,1);                                                
    dcc_vals = zeros(nData,1);                                                

    for i = 1:nData                                                           
        N_curr = N_list(i);                                                   
        M_ss   = M_list(i);                                                    
        O_ss   = O_list(i);                                                   
 
        A = calcJacobian(M_ss, O_ss, N_curr, r, s, dr, dd, a1, a2, a3, a4, a5);

 
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

        fprintf('Row #%d: N=%.2e, tau=%.4g, dCC=%.4g\n', i, N_curr, tau_val, dcc_val);
    end

 
    crit1 = 1*10^14.8999;                       
    crit2 = 1*10^15.6026;                      

 
 
 
    figure;                                                                   
    semilogx(N_list, tau_vals, 'o-','LineWidth',1.2,'MarkerSize',6,'DisplayName','\tau');        
    hold on;                                                                   
    xline(crit1, 'r--','LineWidth',1.5,'DisplayName','低氧临界');            
    xline(crit2, 'g--','LineWidth',1.5,'DisplayName','高氧临界');           
    xlabel('N'); ylabel('\tau');                                               
    title('低稳态：自相关时间 \tau vs N');                                    
    legend('Location','best');                                                
    grid on;                                                                  

 
    figure;                                                                    
    semilogx(N_list, dcc_vals, 's-','LineWidth',1.2,'MarkerSize',6,'DisplayName','\Delta CC');           
    hold on;                                                                   
    xline(crit1, 'r--','LineWidth',1.5,'DisplayName','低氧临界');            
    xline(crit2, 'g--','LineWidth',1.5,'DisplayName','高氧临界');          
    xlabel('N'); ylabel('\Delta CC');                                        
    title('低稳态：\Delta CC vs N');                                       
    legend('Location','best');                                              
    grid on;                                                                
end

 
function A = calcJacobian(M_ss, O_ss, N, r, s, dr, dd, a1, a2, a3, a4, a5)
 
    F1 = @(M,O)(0.5*omega(O,dr,dd).*(N + r) - s.*M - 0.5*phi(O,a1,a2,a3,a4,a5).*M.^0.7); 
    F2 = @(M,O)(    omega(O,dr,dd).*N - (1-omega(O,dr,dd)).*r - s.*M -    phi(O,a1,a2,a3,a4,a5).*M.^0.7); 
    eps_rel = 1e-6;                                                       
    dM = max(M_ss*eps_rel, 1e-8);                                           
    dO = max(O_ss*eps_rel, 1e-8);                                          

 
    dF1_dM = (F1(M_ss+dM,O_ss) - F1(M_ss-dM,O_ss)) / (2*dM);                
    dF2_dM = (F2(M_ss+dM,O_ss) - F2(M_ss-dM,O_ss)) / (2*dM);                 
    dF1_dO = (F1(M_ss,O_ss+dO) - F1(M_ss,O_ss-dO)) / (2*dO);                 
    dF2_dO = (F2(M_ss,O_ss+dO) - F2(M_ss,O_ss-dO)) / (2*dO);              

    A = [dF1_dM, dF1_dO; dF2_dM, dF2_dO];                                    
end

function val = omega(O_val, dr, dd)
 
    gamma = O_val./(dr + O_val);                                           
    delta = O_val./(dd + O_val);                                          
    val   = (1 - gamma).*(1 - delta);                                     
end

function val = phi(O_val, a1,a2,a3,a4,a5)
 
    p   = log10(O_val);                                                    
    val = 10.^(a1*p.^4 + a2*p.^3 + a3*p.^2 + a4*p   + a5);                  
end

function tau_val = calcAutoCorrTime(A)
 
    reVals  = real(eig(A));                                               
    negPart = reVals(reVals < 0);                                        
    if isempty(negPart)                                                      
        tau_val = NaN;                                                      
    else
        tau_val = -1 / max(negPart);                                     
    end
end

function dcc_val = calcDeltaCC(A, Sigma, tMax, nPts)
 
    t_vals = linspace(0, tMax, nPts);                                  
    diffSq = zeros(1, nPts);                              
    for k = 1:nPts                                                        
        Msig      = expm(A * t_vals(k)) * Sigma;                          
        diffSq(k) = (Msig(1,2) - Msig(2,1))^2;                                
    end
    area    = trapz(t_vals, diffSq);                                        
    dcc_val = sqrt(area / tMax);                                         
end
